#' run GWAS
#'
#' @param markers A dataframe with marker information. Row names should be chromosomes, column names
#' genotypes.
#' @param map A dataframe with four columns; chromosome, the number of the chromosome,
#' position, the position of the snp on the chromosome, snp.name, the name of the snp and
#' cum.position, the cumulative position of the snp starting from the first chromosome.
#' @param Y an n x p matrix or dataframe of observed phenotypes, on p traits or environments for n
#' individuals. Missing values are allowed.
#' @param K an n x n kinship matrix.
#' @param X an n x c covariate matrix, c being the number of covariates and n being the number
#' of genotypes. c has to be at least one (typically an intercept). No missing values are allowed.
#' If not provided a vector of 1s is used.
#' @param subsetMarkers should the marker data be subsetted?
#' @param markerSubset numeric or character vector for subsetting the markers. Ignored if
#' subsetMarkers = \code{FALSE}.
#' @param snpCovariates a character vacter of snps that are to be included as covariates.
#' @param MAF a numeric value between 0 and 1. Snps with a minor allele frequency outside MAF
#' and 1 - MAF are excluded from the GWAS analysis.
#' @param fitVarComp should the variance components be fitted? If \code{FALSE} they should be supplied
#' in Vg and Ve
#' @param covModel an integer value for the model used when fitting the variance components.
#' \enumerate{
#' \item{unstructured for both Vg and Ve (as in Zhou and Stephens)}
#' \item{unstructered for both Vg and Ve (pairwise, as in Furlotte and Eskin)}
#' \item{factor-analytic for both Vg and Ve}
#' \item{approximate ...}
#' }
#' Ignored if fitVarComp = \code{FALSE}
#' @param vEDiag Should there be environmental correlations if covModel = 2? If traits are measured on
#' the same individuals put \code{FALSE}.
#' @param tolerance a numerical value. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param maxIter an integer. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param maxDiag a numerical value. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param mG an integer. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param mE an integer. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param CmHet a boolean. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param DmHet a boolean. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param stopIfDecreasing a boolean. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param computeLogLik a boolean. Used when fitting the factor analytical model if covModel = 3.
#' See \code{\link{EMFA}}.
#' @param Vg an optional matrix with genotypic variance components. Vg should have row names
#' column names corresponding to the column names of Y. It may contain additional rows and colums
#' which will be ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param Ve an optional matrix with environmental variance components. Ve should have row names
#' column names corresponding to the column names of Y. It may contain additional rows and colums
#' which will be ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param reduceK if \code{TRUE} the kinship matrix is reduced. See \code{\link{reduceKinship}}
#' @param nPca an integer giving the number of Pcas used whe reducing the kinship matrix.
#' Ignored if reduceK = \code{FALSE}
#'
#' @return a list containing the following items:
#' \itemize{
#' \item{\code{Vg} a matrix with genotypic variance compontents.}
#' \item{\code{Ve} a matrix with environmental variance compontents.}
#' \item{\code{M} a matrix of effect size estimates.}
#' \item{\code{TStat} a matrix of t-statistics.}
#' \item{\code{results} a vector of p-values.}
#' \item{\code{resultsWald} a vector of p-values for the wald test.}
#' \item{\code{MExtended} a matrix of snp information, lod scores, lod scores for the Wald
#' test and effect size estimates.}
#' \item{\code{TStatExtended} a matrix of snp information, lod scores, lod scores for the Wald
#' test and t-statistics.}
#' }
#'
#' @references Dahl et al. (2014). Network inference in matrix-variate Gaussian models with
#' non-indenpent noise.
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear mixed model algorithms for
#' genome-wide association studies. Nature Methods, February 2014, Vol. 11, p. 407–409
#'
#' @export


# #snpCovariates.list <- list('AX-90548584')

## TO DO: more error checking
## TO DO: MAX.DIAG SHOULD DEPEND ON THE SCALE OF THE DATA
## TO DO: the following option is still under construction; leave to zero
## LOD.thr <- 0 if larger than zero, it is assumed a GWAS was done previously with the same name
## .. and GWAS is now only (re)run for markers with -log(p) larger than LOD.thr

runGWAS <- function(markers,
  map,
  Y,
  K,
  X = NULL,
  subsetMarkers = FALSE,
  markerSubset = "",
  snpCovariates = "",
  MAF = 0.05,
  fitVarComp = TRUE,
  covModel = 1,
  vEDiag = TRUE,
  tolerance = 1e-6,
  maxIter = 2e5,
  maxDiag = 1e4,
  mG = 1,
  mE = 1,
  CmHet = TRUE,
  DmHet = TRUE,
  stopIfDecreasing = TRUE,
  computeLogLik = TRUE,
  Vg = NULL,
  Ve = NULL,
  reduceK = FALSE,
  nPca = NULL) {

  if (is.null(markers) || !is.data.frame(markers)) stop("markers should be a dataframe")
  if (is.null(map) || !is.data.frame(map)) stop("map should be a dataframe")
  if (is.null(Y) || !(is.matrix(Y) || is.data.frame(Y))) stop("Y should be a matrix or dataframe")
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (is.null(K) || !is.matrix(K)) stop("K should be a matrix")
  if (!is.null(X) && !(is.matrix(X))) stop("X should either be NULL or a matrix")
  if (anyNA(X)) stop("No missing values allowed in X")
  if (is.null(X)) X <- matrix(data = 1, nrow = nrow(Y), ncol = 1)
  if (subsetMarkers && markerSubset == "") stop("If subsetting markers, markerSubset cannot be empty")
  if(!snpCovariates == "" && !all(snpCovariates %in% rownames(markers)))
    stop("All snpCovariates should be in markers")

  ## Check Vg and Ve if variance components are not fitted.
  if (!fitVarComp) {
    if (is.null(Vg) || !is.matrix(Vg)) stop("Vg should be a matrix")
    if (is.null(Ve) || !is.matrix(Ve)) stop("Ve should be a matrix")
    if (is.null(colnames(Vg)) || is.null(rownames(Vg)) ||
        any(colnames(Vg) != rownames(Vg)) || !all(colnames(Vg) %in% colnames(Y)))
      stop("Column names and rownames of Vg should be identical and included in column names of Y")
    if (is.null(colnames(Ve)) || is.null(rownames(Ve)) ||
        any(colnames(Ve) != rownames(Ve)) || !all(colnames(Ve) %in% colnames(Y)))
      stop("Column names and rownames of Ve should be identical and included in column names of Y")
    Vg <- Vg[colnames(Y), colnames(Y)]
    Ve <- Ve[colnames(Y), colnames(Y)]
    colnames(Vg) <- rownames(Vg) <- NULL
    colnames(Ve) <- rownames(Ve) <- NULL
  } else {
    if (is.null(covModel))
      stop("If variance components are computed, covModel cannot be NULL")
  }

  if (reduceK && is.null(nPca))
    stop("If the kinship matrix is to be reduced, nPca cannot be NULL")

  if (covModel %in% c(2, 4)) {stopifnot(snpCovariates=='')}

  ## Make sure that when subsetting markers snpCovariates are included in the subset
  if (subsetMarkers) {
    if (snpCovariates != "") {
      if (length(which(rownames(markers) %in% snpCovariates)) > 0) {
        markerSubset <- sort(c(markerSubset, which(rownames(markers) %in% snpCovariates)))
        cat('Co-factors have been added to the marker-subset \n')
      }
    }
    ##K <- IBS(markers)[rownames(Y), rownames(Y)]
    K <- K[rownames(Y), rownames(Y)]
    markers <- markers[markerSubset, ]
    map <- map[markerSubset, ]
  }

  ## Add snpCovariates to X
  if (snpCovariates[1] != "") {
    X <- cbind(X, t(as.matrix(markers[snpCovariates, rownames(X)])))
    if (ncol(X) == length(snpCovariates)) {
      XRed <- matrix(nrow = nrow(X), ncol = 0, dimnames = list(rownames(X)))
    } else {
      XRed <- as.matrix(X[, 1:(ncol(X) - length(snpCovariates))])
    }
  }

  if (reduceK) {
    K <- reduceKinship(K = K, nPca = nPca)
  }

  ## fit variance components
  if (fitVarComp) {
    ## Unstructured (pairwise) models
    if (covModel == 2) {
      YLong <- reshape2::melt(data = data.frame(genotype = rownames(Y), Y),
        id.vars = "genotype", variable.name = "trait", value.name = "pheno")

      out <- asreml_unstructured_pairwise(d = YLong, K = K, fix.diag = FALSE,
        correlation.matrix = TRUE,
        vE.diag = vE.diag,
        genotype.column = 1,
        traitname.column = 2,
        phenotype.column = 3,
        covariates = integer())

      out$vG.matrix <- out$vG.matrix[colnames(Y), colnames(Y)]
      out$vE.matrix <- out$vE.matrix[colnames(Y), colnames(Y)]
      out$vG.vector <- out$vG.vector[colnames(Y)]
      out$vE.vector <- out$vE.vector[colnames(Y)]
      vG.matrix <- as.matrix(Matrix::nearPD(out$vG.matrix, corr = TRUE)$mat)
      vE.matrix <- as.matrix(Matrix::nearPD(out$vE.matrix, corr = TRUE)$mat)
      vG.matrix <- (matrix(sqrt(out$vG.vector)) %*% t(matrix(sqrt(out$vG.vector)))) * vG.matrix
      vE.matrix <- (matrix(sqrt(out$vE.vector)) %*% t(matrix(sqrt(out$vE.vector)))) * vE.matrix
      rownames(vG.matrix) <- rownames(vE.matrix) <- colnames(Y)
      colnames(vG.matrix) <- colnames(vE.matrix) <- colnames(Y)
      Vg <- vG.matrix
      Ve <- vE.matrix
    } else if (covModel == 3) {
      ## FA models
      ## Including snpCovariates.
      varcomp <- EMFA(Y = Y,
        K = K,
        X = X,
        maxIter = maxIter,
        tolerance = tolerance,
        mG = mG,
        mE = mE,
        CmHet = CmHet,
        DmHet = DmHet,
        maxDiag = maxDiag,
        stopIfDecreasing = stopIfDecreasing,
        computeLogLik = computeLogLik)
      Vg <- solve(varcomp$Cm)
      Ve <- solve(varcomp$Dm)

      if (snpCovariates[1] != "") {
        ## Without snpCovariates.
        varcompRed <- EMFA(Y = Y,
          K = K,
          X = XRed,
          maxIter = maxIter,
          tolerance = tolerance,
          mG = mG,
          mE = mE,
          CmHet = TRUE,
          DmHet = TRUE,
          maxDiag = maxDiag,
          computeLogLik = computeLogLik,
          stopIfDecreasing = stopIfDecreasing)
        VgRed <- solve(varcompRed$Cm)
        VeRed <- solve(varcompRed$Dm)
      }
    } else if (covModel==4) {
      ## ??
      p <- ncol(Y)
      geno <- rownames(Y)
      GBLUP <- sapply(Y, function(i) {
        outH2 <- heritability::marker_h2_means(data.vector = i, geno.vector = geno, K = K)
        delta <- outH2$va / outH2$ve
        return(delta * K %*% solve(delta * K + diag(p), matrix(i)))})
      Vg <- cov(GBLUP)
      Ve <- cov(Y - GBLUP)
    }
  }

  ## Run GWAS
  w <- eigen(K, symmetric = TRUE)
  Dk <- w$values
  Uk <- w$vectors
  Yt <- t(Y) %*% Uk
  colnames(Yt) <- rownames(Y)
  if (ncol(X) > 0) {
    Xt <- t(X) %*% Uk
  }
  VInvArray <- makeVInvArray(Vg = Vg, Ve = Ve, Dk = Dk)
  if (snpCovariates[1] != "") {
    if (ncol(XRed) > 0) {
      XtRed <- t(XRed) %*% Uk
    }
    VInvArrayRed <- makeVInvArray(Vg = VgRed, Ve = VeRed, Dk = Dk)
  }

  nn <- nrow(map)
  means <- rowMeans(markers[1:nn, rownames(Y)])
  means <- means / max(markers)
  excludedMarkers <- which(means < MAF | means > 1 - MAF)

  if (snpCovariates != "") {
    snpCovariateNumbers <- which(rownames(markers) %in% snpCovariates)
    excludedMarkers <- c(excludedMarkers, snpCovariateNumbers)

    extraExcludedMarkers <- numeric()
    for (snp in snpCovariateNumbers) {
      candidates <- which(means == means[snp])
      ## Only the snp itself is not enough; there needs to be one other snp at least with
      ## the same maf, before proceding.
      if (length(candidates) > 1) {
        snpInfo <- markers[snp, rownames(Y)]
        exclude <- apply(markers[candidates, rownames(Y)], 1,
          function(x) {identical(as.numeric(x), as.numeric(snpInfo))})
        extraExcludedMarkers <- c(extraExcludedMarkers, setdiff(candidates[exclude], snp))
      }
    }
    excludedMarkers <- c(excludedMarkers, extraExcludedMarkers)
    snpCovariateNumbers <- sort(c(snpCovariateNumbers, extraExcludedMarkers))
  }

  ## Scan
  p <- ncol(Y)
  M <- TStat <- matrix(nrow = nn, ncol = p, dimnames = list(rownames(markers), colnames(Y)))
  results <- resultsWald <- rep(NA, nn)
  results[excludedMarkers] <- resultsWald[excludedMarkers] <- 1
  markers <- markers[ , rownames(Y)]

  if (snpCovariates[1]!='') {
    est0Red <- estimateEffects(X = XtRed, Y = Yt, VInvArray = VInvArrayRed, returnAllEffects = TRUE)
    fittedMean0Red <- matrix(est0Red$effects.estimates,ncol = length(est0Red$effects.estimates) / p) %*% XtRed
    SS0Red <- LLQuadFormDiag(Y = Yt - fittedMean0Red, VInvArray = VInvArrayRed)
    for (mrk in snpCovariateNumbers) {
      x <- matrix(as.numeric(markers[mrk, ]))
      xt <- t(x) %*% Uk
      LRTRes <- LRTTest(X = XtRed, x = xt, Y = Yt, VInvArray = VInvArrayRed, SS0 = SS0Red)
      results[mrk] <- LRTRes$pvalue
      resultsWald[mrk] <- pchisq(sum((LRTRes$effects / LRTRes$effects.se) ^ 2),
        df = p, lower.tail = FALSE)
      M[mrk, ] <- LRTRes$effects
      TStat[mrk, ] <-  LRTRes$effects / LRTRes$effects.se
    }
  }

  est0 <- estimateEffects(X = Xt, Y = Yt, VInvArray = VInvArray, returnAllEffects = TRUE)
  fittedMean0 <- matrix(est0$effects.estimates, ncol = length(est0$effects.estimates) / p) %*% Xt
  SS0 <- LLQuadFormDiag(Y = Yt - fittedMean0, VInvArray = VInvArray)
  for (mrk in setdiff(1:nn, excludedMarkers)) {
    x <- matrix(as.numeric(markers[mrk, ]))
    xt <- t(x) %*% Uk
    LRTRes <- LRTTest(X = Xt, x = xt, Y = Yt, VInvArray = VInvArray, SS0 = SS0)
    results[mrk] <- LRTRes$pvalue
    resultsWald[mrk] <- pchisq(q = sum((LRTRes$effects / LRTRes$effects.se) ^ 2),
      df = p, lower.tail = FALSE)
    M[mrk, ] <- LRTRes$effects
    TStat[mrk, ] <- LRTRes$effects / LRTRes$effects.se
    if (mrk %% 500 == 0) {cat("Progress: ", (mrk / nn) * 100, " percent\n")}
  }

  MExtended  <- data.frame(map[c("snp.name", "chromosome", "position")], LOD_F = -log10(results),
    LOD_Wald = -log10(resultsWald), M, row.names = rownames(M))

  TStatExtended  <- data.frame(map[c("snp.name", "chromosome", "position")], LOD_F = -log10(results),
    LOD_Wald = -log10(resultsWald), TStat, row.names = rownames(TStat))
  # }
  return(list(Vg = Vg, Ve = Ve, M = M, TStat = TStat, results = results, resultsWald = resultsWald,
    MExtended = MExtended, TStatExtended = TStatExtended))
}

