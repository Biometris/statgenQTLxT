#' Perform multi-trait GWAS
#'
#' \code{runMultiTraitGwas} performs a multi-trait Genome Wide Association
#' Study (GWAS) on phenotypic and' genotypic data contained in a \code{gData}
#' object.
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param subsetMarkers Should the marker data be subsetted?
#' @param markerSubset A numeric or character vector used for subsetting the
#' markers. Ignored if subsetMarkers = \code{FALSE}.
#' @param fitVarComp Should the variance components be fitted? If \code{FALSE}
#' they should be supplied in Vg and Ve
#' @param covModel An integer value for the model used when fitting the variance
#'  components.
#' \enumerate{
#' \item{unstructured for both Vg and Ve (as in Zhou and Stephens (2014))}
#' \item{unstructered for both Vg and Ve (pairwise, as in Furlotte and Eskin
#' (2013))}
#' \item{factor-analytic for both Vg and Ve}
#' \item{approximate ...(as in Kruijer et al. (2015))}
#' }
#' Ignored if fitVarComp = \code{FALSE}
#' @param VeDiag Should there be environmental correlations if covModel = 1 or
#' 2? If traits are measured on the same individuals put \code{FALSE}.
#' @param tolerance A numerical value. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param maxIter An integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param maxDiag A numerical value. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param mG An integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param mE An integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param CmHet A boolean. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param DmHet A boolean. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param stopIfDecreasing A boolean. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param computeLogLik A boolean. Used when fitting the factor analytical model
#' if covModel = 3. See \code{\link{EMFA}}.
#' @param Vg An optional matrix with genotypic variance components. Vg should
#' have row names column names corresponding to the column names of Y. It may
#' contain additional rows and colums which will be ignored. Ignored if
#' fitVarComp = \code{TRUE}.
#' @param Ve An optional matrix with environmental variance components. Ve
#' should have row names column names corresponding to the column names of Y.
#' It may contain additional rows and colums which will be ignored. Ignored if
#' fitVarComp = \code{TRUE}.
#' @param reduceK Should the kinship matrix be reduced? See
#' \code{\link{reduceKinship}}
#' @param nPca An integer giving the number of Pcas used whe reducing the
#' kinship matrix. Ignored if reduceK = \code{FALSE}.
#' @param estCom Should the common SNP-effect model be fitted?
#' @param parallel Should the computation of variance components be done in
#' parallel? Only used if \code{covModel = 2}. A parallel computing environment
#' has to be setup by the user.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @references Dahl et al. (2013). Network inference in matrix-variate Gaussian
#' models with non-independent noise. arXiv preprint arXiv:1312.1622.
#' @references Furlotte, N.A. and Eskin, E. (2015). Efficient multiple-trait
#' association and estimation of genetic correlation using the matrix-variate
#' linear mixed model. Genetics, May 2015, Vol.200-1, p. 59-68.
#' @references Kruijer et al. (2015) Marker-based estimation of heritability in
#' immortal populations. Genetics. February 2015, Vol. 199-2, p. 379-398.
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409.
#'
#' @export
runMultiTraitGwas <- function(gData,
                              environments = NULL,
                              covar = NULL,
                              snpCov = NULL,
                              kin = NULL,
                              kinshipMethod = c("astle", "GRM", "IBS",
                                                "vanRaden"),
                              GLSMethod = 1,
                              subsetMarkers = FALSE,
                              markerSubset = "",
                              MAF = 0.01,
                              fitVarComp = TRUE,
                              covModel = 1,
                              VeDiag = TRUE,
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
                              nPca = NULL,
                              estCom = FALSE,
                              parallel = FALSE) {
  ## Check input.
  if (missing(gData) || !is.gData(gData) || is.null(gData$markers) ||
      is.null(gData$map) || is.null(gData$pheno)) {
    stop(paste("gData should be a valid gData object containing at least map,",
               "markers and pheno.\n"))
  }
  if (!inherits(gData$markers, "Matrix")) {
    stop(paste("markers in gData should be a numerical matrix. Use",
               "recodeMarkers first for recoding."))
  }
  if (anyNA(gData$markers)) {
    stop("markers contains missing values. Impute or remove these first.\n")
  }
  if (!is.null(environments) && ((!is.numeric(environments) &&
                                  !is.character(environments)) ||
                                 length(environments) > 1)) {
    stop("environments should be a single numeric or character value.\n")
  }
  if ((is.character(environments) &&
       !all(environments %in% names(gData$pheno))) ||
      (is.numeric(environments) && any(environments > length(gData$pheno)))) {
    stop("environments should be list items in pheno.\n")
  }
  if (is.null(environments) && length(gData$pheno) > 1) {
    stop("pheno contains multiple environments. Environment cannot be NULL.\n")
  }
  ## SNPs with MAF == 0 always have to be removed to prevent creation of
  ## singular matrices.
  if (MAF <= 1e-6) {
    MAF <- 1e-6
  }
  ## If environments is null set environments to only environment in pheno.
  if (is.null(environments)) {
    environments <- 1
  }
  markers <- gData$markers
  map <- gData$map
  if (!is.null(covar) && !is.numeric(covar) && !is.character(covar)) {
    stop("covar should be a numeric or character vector.\n")
  }
  if ((is.character(covar) && !all(covar %in% colnames(gData$covar))) ||
      (is.numeric(covar) && any(covar > ncol(gData$covar)))) {
    stop("covar should be columns in covar.\n")
  }
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) {
    covar <- colnames(gData$covar)[covar]
  }
  if (!is.null(snpCov) &&
      !all(snpCov %in% colnames(gData$markers))) {
    stop("All snpCovariates should be in markers.\n")
  }
  if (GLSMethod == 1 && !is.null(kin) && !(inherits(kin, "Matrix") ||
                                           is.matrix(kin))) {
    stop("kin should be a matrix.\n")
  }
  if (GLSMethod == 2 && !is.null(kin) &&
      (!is.list(kin) || !all(sapply(kin, FUN = function(k) {
        is.matrix(k) || inherits(k, "Matrix")})) ||
       length(kin) != length(unique(gData$map$chr)))) {
    stop(paste("kin should be a list of matrices of length equal to the number",
               "of chromosomes in the map.\n"))
  }
  if ((GLSMethod == 1 && is.null(gData$kinship) && is.null(kin)) ||
      GLSMethod == 2 && is.null(kin)) {
    kinshipMethod <- match.arg(kinshipMethod)
  }
  if (subsetMarkers && markerSubset == "") {
    stop("If subsetting markers, markerSubset cannot be empty.\n")
  }
  ## Check Vg and Ve if variance components are not fitted.
  if (!fitVarComp) {
    if (is.null(Vg) || !is.matrix(Vg)) {
      stop("Vg should be a matrix.\n")
    }
    if (is.null(Ve) || !is.matrix(Ve)) {
      stop("Ve should be a matrix.\n")
    }
    if (is.null(colnames(Vg)) || is.null(rownames(Vg)) ||
        any(colnames(Vg) != rownames(Vg)) ||
        !all(colnames(Vg) %in% colnames(gData$pheno[[1]])[-1])) {
      stop(paste("Column names and rownames of Vg should be identical and",
                 "included in column names of Y.\n"))
    }
    if (is.null(colnames(Ve)) || is.null(rownames(Ve)) ||
        any(colnames(Ve) != rownames(Ve)) ||
        !all(colnames(Ve) %in% colnames(gData$pheno[[1]])[-1])) {
      stop(paste("Column names and rownames of Ve should be identical and",
                 "included in column names of pheno.\n"))
    }
    Vg <- Vg[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    Ve <- Ve[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    colnames(Vg) <- rownames(Vg) <- NULL
    colnames(Ve) <- rownames(Ve) <- NULL
  } else {
    if (is.null(covModel)) {
      stop("If variance components are computed, covModel cannot be NULL.\n")
    }
  }
  if (reduceK && is.null(nPca)) {
    stop("If the kinship matrix is to be reduced, nPca cannot be NULL.\n")
  }
  ## Make sure that when subsetting markers snpCovariates are included in
  ## the subset
  if (subsetMarkers) {
    if (!is.null(snpCov)) {
      if (!all(which(colnames(markers) %in% snpCov) %in% markerSubset)) {
        markerSubset <- union(markerSubset,
                              which(colnames(markers) %in% snpCov))
        cat("snpCovariates have been added to the marker-subset\n")
      }
    }
    markersRed <- markers[, markerSubset]
    mapRed <- map[markerSubset, ]
  } else {
    markersRed <- markers[, colnames(markers) %in% rownames(map)]
    mapRed <- map[rownames(map) %in% colnames(markers), ]
  }
  ## Keep option open for extension to multiple environments.
  env <- environments
  ## Add covariates to phenotypic data.
  phExp <- expandPheno(gData = gData, env = env, covar = covar,
                       snpCov = snpCov)
  phEnv <- phExp$phEnv
  covEnv <- phExp$covEnv
  ## Convert pheno and covariates to format suitable for fitting var components.
  X <- cbind(rep(1, nrow(phEnv)), as(as.matrix(phEnv[covEnv]), "dgeMatrix"))
  rownames(X) <- phEnv$genotype
  ## Add snpCovariates to X
  if (!is.null(snpCov)) {
    if (ncol(X) == length(snpCov)) {
      XRed <- Matrix::Matrix(nrow = nrow(X), ncol = 0,
                             dimnames = list(rownames(X)))
    } else {
      XRed <- X[, 1:(ncol(X) - length(snpCov)), drop = FALSE]
    }
  }
  ## Construct Y from pheno data in gData.
  Y <- phEnv[, !colnames(phEnv) %in% covEnv]
  rownames(Y) <- Y[["genotype"]]
  Y <- as(as.matrix(Y[, -which(colnames(Y) == "genotype")]), "dgeMatrix")
  if (anyNA(Y)) {
    stop("Phenotypic data cannot contain any missing values.\n")
  }
  if (GLSMethod == 1) {
    ## Compute kinship matrix.
    K <- computeKin(GLSMethod = 1, kin = kin, gData = gData,
                    markers = markersRed, kinshipMethod = kinshipMethod)
    K <- K[rownames(K) %in% rownames(Y), colnames(K) %in% rownames(Y)]
    if (reduceK) {
      K <- reduceKinship(K = K, nPca = nPca)
    }
    Y <- Y[rownames(Y) %in% rownames(K), ]
    X <- X[rownames(X) %in% rownames(K), , drop = FALSE]
  } else if (GLSMethod == 2) {
    ## Compute kinship matrices per chromosome. Only needs to be done once.
    chrs <- unique(mapRed$chr[rownames(mapRed) %in% colnames(markersRed)])
    KChr <- computeKin(GLSMethod = 2, kin = kin, gData = gData,
                       markers = markersRed, map = mapRed,
                       kinshipMethod = kinshipMethod)
    KChr <- lapply(X = KChr, FUN = function(x) {
      x[rownames(x) %in% rownames(Y), colnames(x) %in% rownames(Y)]
    })
    if (reduceK) {
      KChr <- lapply(X = KChr, FUN = reduceKinship, nPca = nPca)
    }
    Y <- Y[rownames(Y) %in% rownames(KChr[[1]]), ]
    X <- X[rownames(X) %in% rownames(KChr[[1]]), , drop = FALSE]
  }
  ## fit variance components
  if (fitVarComp) {
    if (GLSMethod == 1) {
      if (covModel == 1) {
        ## Unstructured models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- covUnstr(Y = Y, K = K, X = if (ncol(X) == 1) {
          NULL
        } else {
          X[, -1, drop = FALSE]
        }, fixDiag = FALSE, VeDiag = VeDiag)
        if (!is.null(snpCov)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- covUnstr(Y = Y, K = K, X = if (ncol(XRed) == 1) {
            NULL
          } else {
            XRed[, -1, drop = FALSE]
          }, fixDiag = FALSE, VeDiag = VeDiag)
        }
      } else if (covModel == 2) {
        ## Unstructured (pairwise) models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- covPW(Y = Y, K = K, X = if (ncol(X) == 1) {
          NULL
        } else {
          X[, -1, drop = FALSE]
        }, fixDiag = FALSE, corMat = FALSE, parallel = parallel)
        if (!is.null(snpCov)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- covPW(Y = Y, K = K, X = if (ncol(XRed) == 1) {
            NULL
          } else {
            XRed[, -1, drop = FALSE]
          }, fixDiag = FALSE, corMat = FALSE, parallel = parallel)
        }
      } else if (covModel == 3) {
        ## FA models.
        ## Including snpCovariates.
        varComp <- EMFA(Y = Y, K = K, X = X, maxIter = maxIter,
                        tolerance = tolerance, mG = mG, mE = mE, CmHet = CmHet,
                        DmHet = DmHet, maxDiag = maxDiag,
                        stopIfDecreasing = stopIfDecreasing,
                        computeLogLik = computeLogLik)
        if (!is.null(snpCov)) {
          ## Without snpCovariates.
          varCompRed <- EMFA(Y = Y, K = K, X = XRed, maxIter = maxIter,
                             tolerance = tolerance, mG = mG, mE = mE,
                             CmHet = TRUE, DmHet = TRUE, maxDiag = maxDiag,
                             computeLogLik = computeLogLik,
                             stopIfDecreasing = stopIfDecreasing)
        }
      } else if (covModel == 4) {
        ## ??
        geno <- rownames(Y)
        GBLUP <- sapply(as.matrix(Y), function(i) {
          outH2 <- heritability::marker_h2_means(data.vector = i,
                                                 geno.vector = geno,
                                                 K = as.matrix(K))
          delta <- outH2$va / outH2$ve
          return(delta * K %*% solve((delta * K + diag(nrow(Y))), matrix(i)))})
        varComp <- list(Vg = cov(GBLUP), Ve = cov(Y - GBLUP))
      }
      Vg <- varComp$Vg
      Ve <- varComp$Ve
      if (!is.null(snpCov)) {
        VgRed <- varCompRed$Vg
        VeRed <- varCompRed$Ve
      }
    } else if (GLSMethod == 2) {
      if (covModel == 1) {
        ## Unstructured models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          covUnstr(Y = Y, K = KChr[[which(chrs == chr)]],
                   X = if (ncol(X) == 1) NULL else X[, -1, drop = FALSE],
                   fixDiag = FALSE, VeDiag = VeDiag)
        }, simplify = FALSE)
        if (!is.null(snpCov)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            covUnstr(Y = Y, K = KChr[[which(chrs == chr)]],
                     X = if (ncol(XRed) == 1) NULL else
                       XRed[, -1, drop = FALSE], fixDiag = FALSE,
                     VeDiag = VeDiag)
          }, simplify = FALSE)
        }
      } else if (covModel == 2) {
        ## Unstructured (pairwise) models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          covPW(Y = Y, K = KChr[[which(chrs == chr)]],
                X = if (ncol(X) == 1) {
                  NULL
                } else {
                  X[, -1, drop = FALSE]
                }, fixDiag = FALSE, corMat = TRUE, parallel = parallel)
        }, simplify = FALSE)
        if (!is.null(snpCov)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            covPW(Y = Y, K = KChr[[which(chrs == chr)]],
                  X = if (ncol(XRed) == 1) {
                    NULL
                  } else {
                    XRed[, -1, drop = FALSE]
                  }, fixDiag = FALSE, corMat = TRUE, parallel = parallel)
          }, simplify = FALSE)
        }
      } else if (covModel == 3) {
        ## FA models.
        ## Including snpCovariates.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          EMFA(Y = Y, K = KChr[[which(chrs == chr)]], X = X, maxIter = maxIter,
               tolerance = tolerance, mG = mG, mE = mE, CmHet = CmHet,
               DmHet = DmHet, maxDiag = maxDiag,
               stopIfDecreasing = stopIfDecreasing,
               computeLogLik = computeLogLik)
        }, simplify = FALSE)
        if (!is.null(snpCov)) {
          ## Without snpCovariates.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            EMFA(Y = Y, K = KChr[[which(chrs == chr)]], X = XRed,
                 maxIter = maxIter, tolerance = tolerance, mG = mG, mE = mE,
                 CmHet = CmHet, DmHet = DmHet, maxDiag = maxDiag,
                 stopIfDecreasing = stopIfDecreasing,
                 computeLogLik = computeLogLik)
          }, simplify = FALSE)
        }
      } else if (covModel == 4) {
        ## ??
        geno <- rownames(Y)
        GBLUP <- sapply(as.data.frame(Y), function(i) {
          outH2 <- heritability::marker_h2_means(data.vector = i,
                                                 geno.vector = geno, K = K)
          delta <- outH2$va / outH2$ve
          return(delta * K %*% solve((delta * K + diag(nrow(Y))), matrix(i)))})
        varComp <- list(Vg = cov(GBLUP), Ve = cov(Y - GBLUP))
      }
      Vg <- setNames(lapply(X = varComp, FUN = "[[", 1), paste("chr", chrs))
      Ve <- setNames(lapply(X = varComp, FUN = "[[", 2), paste("chr", chrs))
      if (!is.null(snpCov)) {
        VgRed <- lapply(X = varCompRed, FUN = "[[", 1)
        VeRed <- lapply(X = varCompRed, FUN = "[[", 2)
      }
    } #end GLSMethod 2
  } #end varComp
  allFreq <- Matrix::colMeans(markersRed[rownames(Y),
                                         rownames(mapRed)]) / max(markersRed)
  markersRed <- markersRed[rownames(Y), ]
  ## Run GWAS.
  if (GLSMethod == 1) {
    estEffRes <- estEffTot(markers = markersRed, X = X, Y = Y, K = K,
                           XRed = XRed, Vg = Vg, Ve = Ve, snpCov = snpCov,
                           allFreq = allFreq, MAF = MAF, estCom = estCom)
    list2env(estEffRes, envir = environment())
  } else if (GLSMethod == 2) {
    pValues <- pValCom <- pValQtlE <- numeric()
    ## Create an empty matrix with traits as header.
    effs <- effsSe <- effsCom <- effsComSe <- t(gData$pheno[[env]][FALSE, -1])
    for (chr in chrs) {
      mapRedChr <- mapRed[mapRed$chr == chr, ]
      markersRedChr <-
        markersRed[, colnames(markersRed) %in% rownames(mapRedChr),
                   drop = FALSE]
      allFreqChr <- Matrix::colMeans(markersRedChr) / max(markersRedChr)
      snpCovChr <- snpCov[snpCov %in% colnames(markersRedChr)]
      chrNum <- which(chrs == chr)
      estEffRes <- estEffTot(markers = markersRedChr, X = X, Y = Y,
                             K = KChr[[chrNum]], XRed = XRed, Vg = Vg[[chrNum]],
                             Ve = Ve[[chrNum]], snpCov = snpCovChr,
                             allFreq = allFreqChr, MAF = MAF, estCom = estCom)
      pValues <- c(pValues, estEffRes$pValues)
      effs <- cbind(effs, estEffRes$effs)
      effsSe <- cbind(effsSe, estEffRes$effsSe)
      pValCom <- c(pValCom, estEffRes$pValCom)
      effsCom <- c(effsCom, estEffRes$effsCom)
      effsComSe <- c(effsComSe, estEffRes$effsComSe)
      pValQtlE <- c(pValQtlE, estEffRes$pValQtlE)
    }
  }
  ## Convert effs and effsSe to long format and merge.
  effsLong <- reshape2::melt(effs)
  effsSeLong <- reshape2::melt(effsSe)
  effsTot <- merge(effsLong, effsSeLong, by = c("Var1", "Var2"))
  ## Melt creates factors. Reconvert trait and snp to character.
  effsTot$trait <- as.character(effsTot$Var1)
  effsTot$snp <- as.character(effsTot$Var2)
  ## Bind common effects, SE, and pvalues together.
  if (estCom) {
    comDat <- cbind(pValCom, effsCom, effsComSe, pValQtlE)
    comDat <- as.data.frame(comDat)
    comDat$snp <- rownames(comDat)
  }
  ## Set up a data.frame for storing results containing map info and
  ## allele frequencies.
  GWAResult <- data.frame(snp = rownames(mapRed), mapRed, allFreq = allFreq,
                          row.names = rownames(mapRed),
                          stringsAsFactors = FALSE)
  ## Merge the effects and effectsSe to the results.
  GWAResult <- merge(GWAResult, effsTot, by = "snp")
  GWAResult <- merge(GWAResult, data.frame(snp = names(pValues), pValues,
                                           stringsAsFactors = FALSE),
                     by = "snp")
  if (estCom) {
    GWAResult <- merge(GWAResult, comDat, by = "snp")
  }
  GWAResult$LOD <- -log10(GWAResult$pValues)
  ## Select and compute relevant columns.
  relCols <- c("snp", "trait", "chr", "pos", "pValues", "LOD", "value.x",
               "value.y", "allFreq",
               if (estCom) {c("pValCom", "effsCom", "effsComSe",
                              "pValQtlE")})
  GWAResult <- GWAResult[, relCols]
  colnames(GWAResult)[colnames(GWAResult) %in% c("pValues", "value.x",
                                                 "value.y")] <-
    c("pValue", "effect", "effectSe")
  GWAResult <- GWAResult[order(GWAResult$trait, GWAResult$chr, GWAResult$pos), ]
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   MAF = MAF,
                   GLSMethod =
                     factor(GLSMethod, levels = c(1, 2),
                            labels = c("single kinship matrix",
                                       "chromosome specific kinship matrices")),
                   varComp = list(Vg = Vg, Ve = Ve))
  return(createGWAS(GWAResult = setNames(list(GWAResult),
                                         names(gData$pheno)[env]),
                    signSnp = NULL,
                    kin = if (GLSMethod == 1) {
                      K
                    } else {
                      KChr
                    },
                    thr = NULL,
                    GWASInfo = GWASInfo))
}
