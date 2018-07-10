#' Perform multi-trait GWAS
#'
#' \code{runMultiTraitGwas} performs a multi-trait Genome Wide Association
#' Study (GWAS) on phenotypic and' genotypic data contained in a \code{gData}
#' object.
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param subsetMarkers should the marker data be subsetted?
#' @param markerSubset numeric or character vector used for subsetting the
#' markers. Ignored if subsetMarkers = \code{FALSE}.
#' @param fitVarComp should the variance components be fitted? If \code{FALSE}
#' they should be supplied in Vg and Ve
#' @param covModel an integer value for the model used when fitting the variance
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
#' @param tolerance a numerical value. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param maxIter an integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param maxDiag a numerical value. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param mG an integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param mE an integer. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param CmHet a boolean. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param DmHet a boolean. Used when fitting the factor analytical model if
#' covModel = 3. See \code{\link{EMFA}}.
#' @param stopIfDecreasing a boolean. Used when fitting the factor analytical
#' model if covModel = 3. See \code{\link{EMFA}}.
#' @param computeLogLik a boolean. Used when fitting the factor analytical model
#' if covModel = 3. See \code{\link{EMFA}}.
#' @param Vg an optional matrix with genotypic variance components. Vg should
#' have row names column names corresponding to the column names of Y. It may
#' contain additional rows and colums which will be ignored. Ignored if
#' fitVarComp = \code{TRUE}.
#' @param Ve an optional matrix with environmental variance components. Ve
#' should have row names column names corresponding to the column names of Y.
#' It may contain additional rows and colums which will be ignored. Ignored if
#' fitVarComp = \code{TRUE}.
#' @param reduceK if \code{TRUE} the kinship matrix is reduced. See
#' \code{\link{reduceKinship}}
#' @param nPca an integer giving the number of Pcas used whe reducing the
#' kinship matrix. Ignored if reduceK = \code{FALSE}.
#' @param parallel Should the computation of variance components be done in
#' parallel. Only used if \code{covModel = 2}. A parallel computing environment
#' has to be setup by the user.
#'
#' @return an object of class \code{\link{GWAS}}.
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
#' @importFrom rlang .data
#'
#' @export

runMultiTraitGwas <- function(gData,
                              environments = NULL,
                              covar = NULL,
                              snpCovariates = NULL,
                              kin = NULL,
                              kinshipMethod = "astle",
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
  if ((is.character(environments) && !all(environments %in% names(gData$pheno))) ||
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
  ## Construct Y from pheno data in gData. Remove rownames first since
  ## column_to_rownames doesn't overwrite existing rownames.
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
  if (!is.null(snpCovariates) && !all(snpCovariates %in% colnames(gData$markers))) {
    stop("All snpCovariates should be in markers.\n")
  }
  if (GLSMethod == 1 && !is.null(kin) && !(inherits(kin, "Matrix") ||
                                           is.matrix(kin))) {
    stop("kin should be a matrix.\n")
  }
  if (GLSMethod == 2 && !is.null(kin) && (!is.list(kin) ||
                                          !all(sapply(kin, FUN = function(k) {
                                            is.matrix(k) ||
                                              inherits(k, "Matrix")})) ||
                                          length(kin) != dplyr::n_distinct(gData$map$chr))) {
    stop(paste("kin should be a list of matrices of length equal to the number",
               "of chromosomes in the map.\n"))
  }
  if ((GLSMethod == 1 && is.null(gData$kinship) && is.null(kin)) ||
      GLSMethod == 2 && is.null(kin)) {
    kinshipMethod <- match.arg(kinshipMethod,
                               choices = c("astle", "GRM", "IBS", "vanRaden"))
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
  if (covModel %in% c(4)) {
    stopifnot(snpCovariates == "")
  }
  ## Make sure that when subsetting markers snpCovariates are included in the subset
  if (subsetMarkers) {
    if (!is.null(snpCovariates)) {
      if (!all(which(colnames(markers) %in% snpCovariates) %in% markerSubset)) {
        markerSubset <- union(markerSubset,
                              which(colnames(markers) %in% snpCovariates))
        cat('snpCovariates have been added to the marker-subset \n')
      }
    }
    markersRed <- markers[, markerSubset]
    mapRed <- map[markerSubset, ]
  } else {
    markersRed <- markers[, colnames(markers) %in% rownames(map)]
    mapRed <- map[rownames(map) %in% colnames(markers), ]
  }
  ## Keep option open for extension to multiple environments.
  environment <- environments
  ## Add covariates to phenotypic data.
  phenoExp <- expandPheno(gData = gData, environment = environment,
                          covar = covar, snpCovariates = snpCovariates)
  phenoEnvir <- phenoExp$phenoEnvir
  covarEnvir <- phenoExp$covarEnvir
  ## Convert pheno and covariates to format suitable for fitting variance components.
  X <- cbind(rep(1, nrow(phenoEnvir)),
             as(as.matrix(phenoEnvir[covarEnvir]), "dgeMatrix"))
  rownames(X) <- phenoEnvir$genotype
  ## Add snpCovariates to X
  if (!is.null(snpCovariates)) {
    if (ncol(X) == length(snpCovariates)) {
      XRed <- Matrix::Matrix(nrow = nrow(X), ncol = 0,
                             dimnames = list(rownames(X)))
    } else {
      XRed <- X[, 1:(ncol(X) - length(snpCovariates)), drop = FALSE]
    }
  }
  Y <- as(as.matrix(tibble::column_to_rownames(
    tibble::remove_rownames(phenoEnvir[, !colnames(phenoEnvir) %in% covarEnvir]),
    var = "genotype")), "dgeMatrix")
  if (anyNA(Y)) {
    stop("Phenotypic data cannot contain any missing values.\n")
  }
  ## Compute kinship matrix.
  if (GLSMethod == 1) {
    if (is.null(kin)) {
      if (!is.null(gData$kinship)) {
        K <- gData$kinship
      } else {
        K <- do.call(kinshipMethod, list(X = gData$markers))
      }
    } else if (is.matrix(kin)) {
      K <- as(kin, "dsyMatrix")
    } else {
      K <- kin
    }
  } else if (GLSMethod == 2) {
    ## Compute kinship matrices per chromosome. Only needs to be done once.
    chrs <- unique(mapRed$chr[rownames(mapRed) %in% colnames(markersRed)])
    if (!is.null(K)) {
      ## K is supplied. Set KChr to K.
      KChr <- lapply(K, FUN = function(k) {
        if (is.matrix(k)) {
          as(k, "dsyMatrix")
        } else {
          k
        }
      })
    } else {
      ## Compute chromosome specific kinship matrices.
      KChr <- chrSpecKin(gData = createGData(geno = markersRed, map = mapRed),
                         kinshipMethod = kinshipMethod)
    }
  }
  if (GLSMethod == 1) {
    K <- K[rownames(K) %in% rownames(Y), colnames(K) %in% rownames(Y)]
    Y <- Y[rownames(Y) %in% rownames(K), ]
    X <- X[rownames(X) %in% rownames(K), , drop = FALSE]
  } else if (GLSMethod == 2) {
    KChr <- lapply(X = KChr, FUN = function(x) {
      x[rownames(x) %in% rownames(Y), colnames(x) %in% rownames(Y)]
    })
    Y <- Y[rownames(Y) %in% rownames(KChr[[1]]), ]
    X <- X[rownames(X) %in% rownames(KChr[[1]]), , drop = FALSE]
  }
  if (reduceK) {
    K <- reduceKinship(K = K, nPca = nPca)
  }
  ## fit variance components
  if (fitVarComp) {
    if (GLSMethod == 1) {
      if (covModel == 1) {
        ## Unstructured models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- covUnstructured(Y = Y, K = K,
                                   X = if (ncol(X) == 1) {
                                     NULL
                                   } else {
                                     X[, -1, drop = FALSE]
                                   },
                                   fixDiag = FALSE, VeDiag = VeDiag)
        if (!is.null(snpCovariates)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- covUnstructured(Y = Y, K = K,
                                        X = if (ncol(XRed) == 1) {
                                          NULL
                                        } else {
                                          XRed[, -1, drop = FALSE]
                                        },
                                        fixDiag = FALSE, VeDiag = VeDiag)
        }
      } else if (covModel == 2) {
        ## Unstructured (pairwise) models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- covPairwise(Y = Y, K = K,
                               X = if (ncol(X) == 1) NULL else
                                 X[, -1, drop = FALSE],
                               fixDiag = FALSE, corMat = FALSE,
                               parallel = parallel)
        if (!is.null(snpCovariates)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- covPairwise(Y = Y, K = K,
                                    X = if (ncol(XRed) == 1) NULL else
                                      XRed[, -1, drop = FALSE],
                                    fixDiag = FALSE, corMat = FALSE,
                                    parallel = parallel)
        }
      } else if (covModel == 3) {
        ## FA models.
        ## Including snpCovariates.
        varComp <- EMFA(Y = Y, K = K, X = X, maxIter = maxIter,
                        tolerance = tolerance, mG = mG, mE = mE, CmHet = CmHet,
                        DmHet = DmHet, maxDiag = maxDiag,
                        stopIfDecreasing = stopIfDecreasing,
                        computeLogLik = computeLogLik)
        if (!is.null(snpCovariates)) {
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
      if (!is.null(snpCovariates)) {
        VgRed <- varCompRed$Vg
        VeRed <- varCompRed$Ve
      }
    } else if (GLSMethod == 2) {
      if (covModel == 1) {
        ## Unstructured models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          covUnstructured(Y = Y,
                          K = KChr[[chrs == chr]],
                          X = if (ncol(X) == 1) NULL else X[, -1, drop = FALSE],
                          fixDiag = FALSE, VeDiag = VeDiag)
        }, simplify = FALSE)
        if (!is.null(snpCovariates)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            covUnstructured(Y = Y,
                            K = KChr[[chrs == chr]],
                            X = if (ncol(XRed) == 1) NULL else
                              XRed[, -1, drop = FALSE],
                            fixDiag = FALSE, VeDiag = VeDiag)
          }, simplify = FALSE)
        }
      } else if (covModel == 2) {
        ## Unstructured (pairwise) models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          covPairwise(Y = Y,
                      K = KChr[[chrs == chr]],
                      X = if (ncol(X) == 1) NULL else X[, -1, drop = FALSE],
                      fixDiag = FALSE, corMat = TRUE, parallel = parallel)
        }, simplify = FALSE)
        if (!is.null(snpCovariates)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            covPairwise(Y = Y,
                        K = KChr[[chrs == chr]],
                        X = if (ncol(XRed) == 1) NULL else
                          XRed[, -1, drop = FALSE],
                        fixDiag = FALSE, corMat = TRUE, parallel = parallel)
          }, simplify = FALSE)
        }
      } else if (covModel == 3) {
        ## FA models.
        ## Including snpCovariates.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          EMFA(Y = Y, K = KChr[[chrs == chr]], X = X, maxIter = maxIter,
               tolerance = tolerance, mG = mG, mE = mE, CmHet = CmHet,
               DmHet = DmHet, maxDiag = maxDiag,
               stopIfDecreasing = stopIfDecreasing,
               computeLogLik = computeLogLik)
        }, simplify = FALSE)
        if (!is.null(snpCovariates)) {
          ## Without snpCovariates.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            EMFA(Y = Y, K = KChr[[chrs == chr]], X = XRed, maxIter = maxIter,
                 tolerance = tolerance, mG = mG, mE = mE, CmHet = CmHet,
                 DmHet = DmHet, maxDiag = maxDiag,
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
      Vg <- setNames(lapply(X = varComp, FUN = function(x) {x[[1]]}),
                     paste("chr", chrs))
      Ve <- setNames(lapply(X = varComp, FUN = function(x) {x[[2]]}),
                     paste("chr", chrs))
      if (!is.null(snpCovariates)) {
        VgRed <- lapply(X = varCompRed, FUN = function(x) {x[[1]]})
        VeRed <- lapply(X = varCompRed, FUN = function(x) {x[[2]]})
      }
    }
  }
  ## Create data.frame and matrices for storing GWAS Results.
  nn <- nrow(mapRed)
  allFreq <- Matrix::colMeans(markersRed[rownames(Y),
                                         rownames(mapRed)]) / max(markersRed)
  effects <- effectsSe <- matrix(nrow = nn, ncol = ncol(Y),
                                 dimnames = list(colnames(markersRed),
                                                 colnames(Y)))
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
                          row.names = rownames(mapRed),
                          stringsAsFactors = FALSE)
  ## Run GWAS.
  if (GLSMethod == 1) {
    segMarkers <- which(allFreq < MAF | allFreq > 1 - MAF)
    ## Add snpCovariates to segregating markers.
    excludedMarkers <- union(c(segMarkers, ncol(markersRed) + 1),
                             computeExcludedMarkers(snpCovariates = snpCovariates,
                                                    markersRed = markersRed,
                                                    allFreq = allFreq))
    if (!is.null(snpCovariates)) {
      effEstSnpCov <- estimateEffects(Y = Y, W = XRed,
                                      X = markersRed[, snpCovariates,
                                                     drop = FALSE],
                                      Vg = Vg, Ve = Ve, K = K)
    }
    effEst <- estimateEffects(Y = Y, W = X,
                              X = markersRed[, -excludedMarkers],
                              Vg = Vg, Ve = Ve, K = K)
    pValues <- c(effEst$pVals,
                 if (!is.null(snpCovariates)) effEstSnpCov$pVals)
    effects <- cbind(effEst$effects,
                     if (!is.null(snpCovariates)) effEstSnpCov$effects)
    effectsSe <- cbind(effEst$effectsSe,
                       if (!is.null(snpCovariates)) effEstSnpCov$effectsSe)
  } else if (GLSMethod == 2) {
    pValues <- numeric()
    ## Create an empty matrix with traits as header.
    effects <- effectsSe <- t(gData$pheno[[environment]][FALSE, -1])
    for (chr in chrs) {
      w <- eigen(KChr[[chrs == chr]], symmetric = TRUE)
      Dk <- w$values
      Uk <- w$vectors
      Yt <- Matrix::crossprod(Y, Uk)
      colnames(Yt) <- rownames(Y)
      if (ncol(X) > 0) {
        Xt <- Matrix::crossprod(X, Uk)
      }
      VInvArray <- makeVInvArray(Vg = Vg[[chrs == chr]], Ve = Ve[[chrs == chr]],
                                 Dk = Dk)
      if (!is.null(snpCovariates)) {
        if (ncol(XRed) > 0) {
          XtRed <- Matrix::crossprod(XRed, Uk)
        }
        VInvArrayRed <- makeVInvArray(Vg = VgRed[[chrs == chr]],
                                      Ve = VeRed[[chrs == chr]], Dk = Dk)
      }
      mapRedChr <- mapRed[mapRed$chr == chr, ]
      markersRedChr <- markersRed[, colnames(markersRed) %in% rownames(mapRedChr),
                                  drop = FALSE]
      allFreqChr <- Matrix::colMeans(markersRedChr) / max(markersRedChr)
      segMarkers <- which(allFreqChr < MAF | allFreqChr > 1 - MAF)
      snpCovChr <- snpCovariates[snpCovariates %in% colnames(markersRedChr)]
      ## Add snpCovariates to segregating markers.
      excludedMarkers <- union(c(segMarkers, ncol(markersRed) + 1),
                               computeExcludedMarkers(snpCovariates = snpCovChr,
                                                      markersRed = markersRedChr,
                                                      allFreq = allFreqChr))
      if (length(snpCovChr) > 0) {
        effEstSnpCov <- estimateEffects(Y = Y, W = XRed,
                                        X = markersRedChr[, snpCovChr,
                                                          drop = FALSE],
                                        Vg = Vg[[chrs == chr]],
                                        Ve = Ve[[chrs == chr]],
                                        K = KChr[[chrs == chr]])
      }
      effEst <- estimateEffects(Y = Y, W = X,
                                X = markersRedChr[, -excludedMarkers],
                                Vg = Vg[[chrs == chr]], Ve = Ve[[chrs == chr]],
                                K = KChr[[chrs == chr]])
      pValues <- c(pValues, effEst$pVals,
                   if (length(snpCovChr) > 0) effEstSnpCov$pVals)
      effects <- cbind(effects, effEst$effects,
                       if (length(snpCovChr) > 0) effEstSnpCov$effects)
      effectsSe <- cbind(effectsSe, effEst$effectsSe,
                         if (length(snpCovChr) > 0) effEstSnpCov$effectsSe)
    }
  }
  ## Convert effects and effectsSe to long format and merge.
  effectsTot <- reshape2::melt(effects) %>%
    dplyr::inner_join(reshape2::melt(effectsSe), by = c("Var1", "Var2")) %>%
    dplyr::mutate(Var1 = as.character(.data$Var1),
                  Var2 = as.character(.data$Var2))
  ## Merge the effects and effectsSe to the results
  GWAResult <- dplyr::inner_join(GWAResult, effectsTot,
                                 by = c("snp" = "Var2")) %>%
    dplyr::inner_join(data.frame(snp = names(pValues), pValues,
                                 stringsAsFactors = FALSE), by = "snp") %>%
    dplyr::mutate(LOD = -log10(.data$pValues)) %>%
    ## Select and compute relevant columns.
    ## Melt creates factors. Reconvert trait to character.
    dplyr::select(.data$snp, trait = .data$Var1, .data$chr, .data$pos,
                  pValue = .data$pValues, .data$LOD, effect = .data$value.x,
                  effectSe = .data$value.y, .data$allFreq) %>%
    dplyr::arrange(.data$trait, .data$chr, .data$pos)
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   MAF = MAF,
                   GLSMethod = factor(GLSMethod, levels = c(1, 2),
                                      labels = c("single kinship matrix",
                                                 "chromosome specific kinship matrices")),
                   varComp = list(Vg = Vg, Ve = Ve))
  return(createGWAS(GWAResult = setNames(list(GWAResult),
                                         names(gData$pheno)[environment]),
                    signSnp = NULL,
                    kin = if (GLSMethod == 1) {
                      if (is.null(kin)) {
                        gData$kinship
                      } else {
                        K
                      }
                    } else {
                      KChr
                    },
                    thr = NULL,
                    GWASInfo = GWASInfo))
}

