#' run multi trait GWAS
#'
#' @param gData an object of class \code{gData} containing at least \code{map}, \code{markers} and
#' \code{pheno}.
#' @param environments a vector of environments on which to run GWAS. These can be either numeric indices or
#' character names of list items in \code{pheno}. If \code{NULL} GWAS is run for all environments.
#' @param K an optional kinship matrix. If \code{NULL} then matrix \code{kinship} in \code{gData}
#' is used. If both \code{K} is provided and \code{gData} contains a matrix \code{kinship}
#' then \code{K} is used.
#' @param covar an optional vector of covariates taken into account when running GWAS. These can be either
#' numeric indices or character names of columns in \code{covar} in \code{gData}. If \code{NULL} no
#' covariates are used.
#' @param snpCovariates an optional character vector of snps to be included as covariates.
#' @param subsetMarkers should the marker data be subsetted?
#' @param markerSubset numeric or character vector for subsetting the markers. Ignored if
#' subsetMarkers = \code{FALSE}.
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

runMultiTraitGwas <- function(gData,
  environments = NULL,
  covar = NULL,
  snpCovariates = NULL,
  K = NULL,
  subsetMarkers = FALSE,
  markerSubset = "",
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
  ## Check input.
  if (missing(gData) || !is.gData(gData) || is.null(gData$markers) ||
      is.null(gData$map) || is.null(gData$pheno))
    stop("gData should be a valid gData object containing at least map, markers and pheno.\n")
  if (!is.null(environments) && !is.numeric(environments) && !is.character(environments))
    stop("environments should be a numeric or character vector.\n")
  if ((is.character(environments) && !all(environments %in% names(gData$pheno))) ||
      (is.numeric(environments) && any(environments > length(gData$pheno))))
    stop("environments should be list items in pheno.\n")
  ## If environments is null set environments to all environments in pheno.
  if (is.null(environments)) environments <- 1:length(gData$pheno)
  markers <- gData$markers
  map <- gData$map
  ## Construct Y from pheno data in gData. Remove rownames first since column_to_rownames doesn't
  ## overwrite existing rownames.
  if (!is.null(covar) && !is.numeric(covar) && !is.character(covar))
    stop("covar should be a numeric or character vector.\n")
  if ((is.character(covar) && !all(covar %in% colnames(gData$covar))) ||
      (is.numeric(covar) && any(covar > ncol(gData$covar))))
    stop("covar should be columns in covar.\n")
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) covar <- colnames(gData$covar)[covar]
  if (!is.null(snpCovariates) && !all(snpCovariates %in% colnames(gData$markers)))
    stop("All snpCovariates should be in markers.\n")
  if (!is.null(K) && !is.matrix(K))
    stop("K should be a matrix")
  if (is.null(K) && is.null(gData$kinship))
    stop("gData contains no matrix kinship so K should be provided.\n")
  if (subsetMarkers && markerSubset == "")
    stop("If subsetting markers, markerSubset cannot be empty.\n")
  ## Check Vg and Ve if variance components are not fitted.
  if (!fitVarComp) {
    if (is.null(Vg) || !is.matrix(Vg))
      stop("Vg should be a matrix.\n")
    if (is.null(Ve) || !is.matrix(Ve))
      stop("Ve should be a matrix.\n")
    if (is.null(colnames(Vg)) || is.null(rownames(Vg)) ||
        any(colnames(Vg) != rownames(Vg)) || !all(colnames(Vg) %in% colnames(gData$pheno[[1]])[-1]))
      stop("Column names and rownames of Vg should be identical and included in column names of Y.\n")
    if (is.null(colnames(Ve)) || is.null(rownames(Ve)) ||
        any(colnames(Ve) != rownames(Ve)) || !all(colnames(Ve) %in% colnames(gData$pheno[[1]])[-1]))
      stop("Column names and rownames of Ve should be identical and included in column names of pheno.\n")
    Vg <- Vg[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    Ve <- Ve[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    colnames(Vg) <- rownames(Vg) <- NULL
    colnames(Ve) <- rownames(Ve) <- NULL
  } else {
    if (is.null(covModel))
      stop("If variance components are computed, covModel cannot be NULL.\n")
  }
  if (reduceK && is.null(nPca))
    stop("If the kinship matrix is to be reduced, nPca cannot be NULL.\n")
  if (is.null(K)) K <- gData$kinship
  if (covModel %in% c(2, 4)) {stopifnot(snpCovariates == "")}
  ## Make sure that when subsetting markers snpCovariates are included in the subset
  if (subsetMarkers) {
    if (!is.null(snpCovariates)) {
      if (!all(length(which(colnames(markers) %in% snpCovariates) %in% markerSubset))) {
        markerSubset <- union(markerSubset,
          which(colnames(markers) %in% snpCovariates))
        cat('snpCovariates have been added to the marker-subset \n')
      }
    }
    markersRed <- markers[, markerSubset]
    mapRed <- map[markerSubset, ]
  } else {
    markersRed <- markers
    mapRed <- map
  }
  environment <- 1
  if (is.null(covar)) {
    phenoEnvir <- gData$pheno[[environment]]
    covarEnvir <- NULL
  } else {
    ## Append covariates to pheno data. Merge to remove values from pheno that are missing in covar.
    phenoEnvir <- merge(gData$pheno[[environment]], gData$covar[covar], by.x = "genotype", by.y = "row.names")
    ## Expand covariates that are a factor (i.e. dummy variables are created) using model.matrix
    ## The new dummies are attached to phenoEnvir, and covar is changed accordingly
    factorCovs <- which(sapply(gData$covar[covar], FUN = is.factor))
    if (length(factorCovs) > 0) {
      covFormula <- as.formula(paste("genotype ~ ", paste(covar[factorCovs], collapse = "+")))
      ## Create dummy variables. Remove intercept.
      extraCov <- as.data.frame(suppressWarnings(model.matrix(object = covFormula, data = phenoEnvir))[, -1])
      ## Add dummy variables to pheno data.
      phenoEnvir <- cbind(phenoEnvir[, -which(colnames(phenoEnvir) %in% names(factorCovs))], extraCov)
      ## Modify covar to suit newly defined columns
      covarEnvir <- c(covar[-factorCovs], colnames(extraCov))
    } else {
      covarEnvir <- covar
    }
  }
  if (!is.null(snpCovariates)) {
    ## Add snp covariates to covar.
    covarEnvir <- c(covarEnvir, snpCovariates)
    ## Add snp covariates to pheno data.
    phenoEnvir <- merge(phenoEnvir, gData$markers[, snpCovariates], by.x = "genotype",
      by.y = "row.names")
    colnames(phenoEnvir)[(ncol(phenoEnvir) - length(snpCovariates) + 1):ncol(phenoEnvir)] <- snpCovariates
  }
  X <- cbind(rep(1, nrow(phenoEnvir)), as.matrix(phenoEnvir[covarEnvir]))
  rownames(X) <- phenoEnvir$genotype
  Y <- as.matrix(tibble::column_to_rownames(tibble::remove_rownames(gData$pheno[[1]]),
    var = "genotype"))
  K <- K[rownames(Y), rownames(Y)]
  ## Add snpCovariates to X
  if (!is.null(snpCovariates)) {
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
      varcomp <- covPairwise(Y = Y, K = K, fixDiag = FALSE, corMat = TRUE, VeDiag = FALSE)
      Vg <- varcomp$Vg
      Ve <- varcomp$Ve
      if (!is.null(snpCovariates)) {
        varcompRed <- covPairwise(Y = Y, K = K, X = X, fixDiag = FALSE, corMat = TRUE, VeDiag = FALSE)
        VgRed <- varcomp$Vg
        VeRed <- varcomp$Ve
      }
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
      colnames(Vg) <- rownames(Vg) <- colnames(Y)
      colnames(Ve) <- rownames(Ve) <- colnames(Y)
      if (!is.null(snpCovariates)) {
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
      p <- nrow(Y)
      geno <- rownames(Y)
      GBLUP <- sapply(as.data.frame(Y), function(i) {
        outH2 <- heritability::marker_h2_means(data.vector = i, geno.vector = geno, K = K)
        delta <- outH2$va / outH2$ve
        return(delta * K %*% solve((delta * K + diag(p)), matrix(i)))})
      Vg <- cov(GBLUP)
      Ve <- cov(Y - GBLUP)
    }
  }
  ## Run GWAS
  w <- eigen(K, symmetric = TRUE)
  Dk <- w$values
  Uk <- w$vectors
  Yt <- crossprod(Y, Uk)
  colnames(Yt) <- rownames(Y)
  if (ncol(X) > 0) {
    Xt <- crossprod(X, Uk)
  }
  VInvArray <- makeVInvArray(Vg = Vg, Ve = Ve, Dk = Dk)
  if (!is.null(snpCovariates)) {
    if (ncol(XRed) > 0) {
      XtRed <- crossprod(XRed, Uk)
    }
    VInvArrayRed <- makeVInvArray(Vg = VgRed, Ve = VeRed, Dk = Dk)
  }
  nn <- nrow(mapRed)
  allFreq <- colMeans(markersRed[rownames(Y), 1:nn]) / max(markersRed)
  excludedMarkers <- which(allFreq < MAF | allFreq > 1 - MAF)
  if (!is.null(snpCovariates)) {
    snpCovariateNumbers <- which(colnames(markersRed) %in% snpCovariates)
    excludedMarkers <- c(excludedMarkers, snpCovariateNumbers)
    extraExcludedMarkers <- numeric()
    for (snp in snpCovariateNumbers) {
      candidates <- which(allFreq == allFreq[snp])
      ## Only the snp itself is not enough; there needs to be one other snp at least with
      ## the same maf, before proceding.
      if (length(candidates) > 1) {
        snpInfo <- markersRed[rownames(Y), snp]
        exclude <- apply(markersRed[rownames(Y), candidates], 2,
          function(x) {identical(as.numeric(x), as.numeric(snpInfo))})
        extraExcludedMarkers <- c(extraExcludedMarkers, setdiff(candidates[exclude], snp))
      }
    }
    excludedMarkers <- c(excludedMarkers, extraExcludedMarkers)
    snpCovariateNumbers <- sort(c(snpCovariateNumbers, extraExcludedMarkers))
  }
  ## Scan
  p <- ncol(Y)
  effects <- effectsSe <- matrix(nrow = nn, ncol = p, dimnames = list(colnames(markersRed), colnames(Y)))
  markersRed <- markersRed[rownames(Y), ]
  GWAResult <- data.frame(trait = NA,
    snp = rownames(mapRed),
    mapRed,
    pValue = NA,
    LOD = NA,
    effect = NA,
    effectSe = NA,
    pValueWald = NA,
    LODWald = NA,
    allFreq = allFreq,
    stringsAsFactors = FALSE)
  if (!is.null(snpCovariates)) {
    est0Red <- estimateEffects(X = XtRed, Y = Yt, VInvArray = VInvArrayRed, returnAllEffects = TRUE)
    fittedMean0Red <- matrix(est0Red$effects.estimates,ncol = length(est0Red$effects.estimates) / p) %*% XtRed
    SS0Red <- LLQuadFormDiag(Y = Yt - fittedMean0Red, VInvArray = VInvArrayRed)
    for (mrk in snpCovariateNumbers) {
      x <- matrix(as.numeric(markersRed[, mrk]))
      xt <- crossprod(x, Uk)
      LRTRes <- LRTTest(X = XtRed, x = xt, Y = Yt, VInvArray = VInvArrayRed, SS0 = SS0Red)
      GWAResult[mrk, "pValue"] <- LRTRes$pvalue
      GWAResult[mrk, "pValueWald"] <- pchisq(sum((LRTRes$effects / LRTRes$effects.se) ^ 2),
        df = p, lower.tail = FALSE)
      effects[mrk, ] <- LRTRes$effects
      effectsSe[mrk, ] <- LRTRes$effects.se
    }
  }
  est0 <- estimateEffects(X = Xt, Y = Yt, VInvArray = VInvArray, returnAllEffects = TRUE)
  fittedMean0 <- matrix(est0$effects.estimates, ncol = length(est0$effects.estimates) / p) %*% Xt
  SS0 <- LLQuadFormDiag(Y = Yt - fittedMean0, VInvArray = VInvArray)
  for (mrk in setdiff(1:nn, excludedMarkers)) {
    x <- matrix(as.numeric(markersRed[, mrk]))
    xt <- crossprod(x, Uk)
    LRTRes <- LRTTest(X = Xt, x = xt, Y = Yt, VInvArray = VInvArray, SS0 = SS0)
    GWAResult[mrk, "pValue"] <- LRTRes$pvalue
    GWAResult[mrk, "pValueWald"] <- pchisq(sum((LRTRes$effects / LRTRes$effects.se) ^ 2),
      df = p, lower.tail = FALSE)
    effects[mrk, ] <- LRTRes$effects
    effectsSe[mrk, ] <-  LRTRes$effects.se
    if (mrk %% 1000 == 0) {cat("Progress: ", (mrk / nn) * 100, " percent\n")}
  }
  ## Add LOD-scores to result
  GWAResult$LOD <- -log10(GWAResult$pValue)
  GWAResult$LODWald <- -log10(GWAResult$pValueWald)
  ## Convert effects en effectsSe to long-format and merge.
  effectsLong <- reshape2::melt(effects)
  effectsSeLong <- reshape2::melt(effectsSe)
  effectsTot <- merge(effectsLong, effectsSeLong, by = c("Var1", "Var2"))
  ## Merge the effects and effectsSe to the results
  GWAResult <- merge(GWAResult, effectsTot, by.x = "snp", by.y = "Var1", sort = FALSE)
  ## Melt creates factors. Reconvert trait to character
  GWAResult$trait <- as.character(GWAResult$Var2)
  GWAResult[c("effect", "effectSe")] <-  GWAResult[c("value.x", "value.y")]
  ## Sort by trait, chr and pos and drop unneeded columns.
  GWAResult <- GWAResult[order(GWAResult$trait, GWAResult$chr, GWAResult$pos),
    !colnames(GWAResult) %in% c("Var2", "value.x", "value.y")]
  ## Remove rownames.
  rownames(GWAResult) <- NULL
  ## Collect info.
  GWASInfo <- list(call = match.call(),
    MAF = MAF,
    varComp = list(Vg = Vg, Ve = Ve))
  return(createGWAS(GWAResult = GWAResult, signSnp = NULL, thr = NULL, GWASInfo = GWASInfo))
}

