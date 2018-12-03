#' Perform multi-trait GWAS
#'
#' \code{runMultiTraitGwas} performs multi-trait or multi-environment Genome
#' Wide Association mapping on phenotypic and genotypic data contained in a
#' \code{gData} object.
#'
#' @section Details: runMultiTraitGwas estimates the effect of a SNP in
#' different environments or on different traits, for each SNP in turn. Genetic
#' and residual covariances are fitted only once, for a model without SNPs.
#' Following the diagonalization scheme of Zhou and Stephens (2014), the
#' following model is fit
#'
#' \eqn{Y = \left(\begin{array}{c} Y_1 \\ \vdots \\ Y_p\end{array}\right) =
#' \left(\begin{array}{c} X_1\gamma_1 \\ \vdots \\ X_p\gamma_p\end{array}\right)
#' + \left(\begin{array}{c} x_1\beta_1 \\ \vdots \\ x_p\beta_p\end{array}\right)
#' + \left(\begin{array}{c} G_1 \\ \vdots \\ G_p\end{array}\right) +
#' \left(\begin{array}{c} E_1 \\ \vdots \\ E_p\end{array}\right)}
#'
#' where Y is a np x 1 vector of phenotypic values for n genotypes and p traits
#' or environments. x is the n x 1 vector of scores for the marker under
#' consideration, and X the n x q design matrix for the other covariates. By
#' default only a trait (environment) specific intercept is included. The vector
#' of genetic background effects
#' (\eqn{\left(\begin{array}{c}G_1 \\ \vdots \\ G_p\end{array}\right)}) is
#' Gaussian with
#' zero mean and covariance \eqn{V_g \otimes K}, where \eqn{V_g} is a p x p
#' matrix of genetic (co)variances, and K an n x n kinship matrix. Similarly,
#' the residual errors
#' (\eqn{\left(\begin{array}{c}E_1 \\ \vdots \\ E_p\end{array}\right)})
#' have covariance
#' \eqn{V_e \otimes I_n}, for a p x p matrix \eqn{V_e} of residual
#' (co)variances.
#'
#' @section Hypotheses for the SNP-effects:
#' For each SNP the null-hypothesis \eqn{\beta_1 = \dots = \beta_p = 0} is
#' tested, using the likelihood ratio test (LRT) described in Zhou and
#' Stephens (2014). If estCom = TRUE, additional tests for a common effect and
#' for QTL x E are performed, using the parameterization \eqn{\beta_j = \alpha +
#' \alpha_j (1 \leq j \leq p)}. As in Korte et al (2012) we use likelihood ratio
#' tests, but not restricted to the bivariate case. For the common effect, we
#' fit the reduced model \eqn{\beta_j = \alpha}, and test if \eqn{\alpha = 0}.
#' For QTL by environment interaction, we test if \eqn{\alpha_1 = \dots =
#' \alpha_p = 0}.
#'
#' @section Models for the genetic and residual covariance:
#' \eqn{V_g} and \eqn{V_e} can be provided by the user
#' (\code{fitVarComp = FALSE});
#' otherwise one of the following models is used, depending on covModel.
#' If \code{covModel = "unst"} an unstructured model is assumed, as in Zhou and
#' Stephens (2014): \eqn{V_g} and \eqn{V_e} can be any positive-definite matrix,
#' requiring a total of p(p+1)/2 parameters per matrix.
#' If \code{covModel = "fa"}, a factor-analytic model is fitted using an
#' EM-algorithm, as in Millet et al (2016). \eqn{V_g} and \eqn{V_e} are assumed
#' to be of the form \eqn{W W^t + D}, where W is a p x m matrix of factor
#' loadings and D a diagonal matrix with trait or environment specific values.
#' m is the order of the model, and the parameters \code{mG} and \code{mE}
#' specify the order used for respectively \eqn{V_g} and \eqn{V_e}.
#' \code{maxIter} sets the maximum number of iterations used in the
#' EM-algorithm.
#' Finally, if \code{covModel = "pw"}, \eqn{V_g} and \eqn{V_e} are estimated
#' 'pairwise', as in Furlotte and Eskin (2015). Looping over pairs of traits
#' or environments \eqn{1 \leq j < k \leq p},
#' \eqn{V_g[j,k] = V_g[k,j]} and \eqn{V_e[j,k] = V_e[k,j]}
#' are estimated assuming a bivariate mixed model. The diagonals of
#' \eqn{V_g} and \eqn{V_e} are fitted assuming univariate mixed models. If the
#' resulting \eqn{V_g} or \eqn{V_e} is not positive-definite, they are
#' replaced by the nearest positive-definite matrix.
#' In case \code{covModel = "unst"} or \code{"pw"} it is possible to assume that
#' \eqn{V_e} is diagonal (\code{VeDiag = TRUE})
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}. The latter should not contain missing
#' values. Multi-trait or multi-environment GWAS is performed for all variables
#' in \code{pheno}.
#' @param environments A vector specifying the environment on which to run GWAS.
#' Thise can be either a numeric index or a character name of a list item in
#' \code{pheno}.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL} no covariates are
#' used. An intercept is included automatically (and should not be assigned as
#' covariate). SNP-covariates should be assigned using the snpCov parameter.
#' @param snpCov An optional character vector of SNP-names to be included as
#' covariates. SNP-names should match those used in \code{gData}.
#' @param subsetMarkers Should GWAS be performed on a subset of markers?
#' @param markerSubset A numeric or character vectorof marker names, indicating
#' on which markers in \code{gData$markers} GWAS is to be performed. Ignored if
#' \code{subsetMarkers = FALSE}.
#' @param MAF The minor allele frequency (MAF) threshold used in GWAS. A
#' numerical value between 0 and 1. SNPs with MAF below this value are not taken
#' into account in the analysis, i.e. p-values and effect sizes are put to
#' missing (\code{NA}).
#' @param fitVarComp Should the variance components be fitted? If \code{FALSE}
#' they should be supplied in Vg and Ve
#' @param covModel A character string indicating the covariance model for the
#' genetic background (Vg) and residual effects (Ve); see details.
#' Either \code{unst} for unstructured for both Vg and
#' Ve (as in Zhou and Stephens (2014)), \code{pw} for unstructered for both Vg
#' and Ve (pairwise, as in Furlotte and Eskin (2013)) or \code{fa} for
#' factor-analytic for both Vg and Ve.\cr
#' Ignored if \code{fitVarComp} = \code{FALSE}
#' @param VeDiag Should there be environmental correlations if covModel = "unst"
#' or "pw"? If traits are measured on the same individuals put \code{TRUE}.
#' @param tolerance A numerical value. The iterating process stops if the
#' difference in conditional log-likelihood between two consecutive iterations
#' drops below \code{tolerance}. Only used when \code{covModel = "fa"}.
#' @param maxIter An integer for the maximum number of iterations. Only used
#' when \code{covModel = "fa"}.
#' @param maxDiag A numical value. The maximal value of the diagonal elements
#' in the precision matrices Cm and Dm (ignoring the low-rank part W W^t).
#' Only used when \code{covModel = "fa"}.
#' @param mG An integer. The order of the genetic part of the factor analytic
#' model. Only used when \code{covModel = "fa"}.
#' @param mE An integer. The order of the environmental part of the factor
#' analytic model. Only used when \code{covModel = "fa"}.
#' @param CmHet Should an extra diagonal part be added in the model for the
#' precision matrix Cm? Only used when \code{covModel = "fa"}.
#' @param DmHet Should an extra diagonal part be added in the model for the
#' precision matrix Dm? Only used when \code{covModel = "fa"}.
#' @param stopIfDecreasing Should the iterating process in the factor analytic
#' model stop if after 50 iterations the log-likelihood decreases between two
#' consecutive iterations? Only used when \code{covModel = "fa"}.
#' @param Vg An optional matrix with genotypic variance components. Vg should
#' have row names column names corresponding to the column names of
#' \code{gData$pheno}. It may contain additional rows and columns which will be
#' ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param Ve An optional matrix with environmental variance components. Ve
#' should have row names column names corresponding to the column names of
#' \code{gData$pheno}. It may contain additional rows and columns which will be
#' ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param reduceK Should the kinship matrix be reduced? See
#' \code{\link{reduceKinship}}
#' @param nPca An integer giving the number of Pcas used when reducing the
#' kinship matrix. Ignored if \code{reduceK} = \code{FALSE}.
#' @param estCom Should the common SNP-effect model be fitted? If \code{TRUE}
#' not only the SNP-effects but also the common SNP-effect and QTL x E effect
#' are estimated.
#' @param parallel Should the computation of variance components be done in
#' parallel? Only used if \code{covModel = "pw"}. A parallel computing
#' environment has to be setup by the user.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @references Dahl et al. (2013). Network inference in matrix-variate Gaussian
#' models with non-independent noise. arXiv preprint arXiv:1312.1622.
#' @references Furlotte, N.A. and Eskin, E. (2015). Efficient multiple-trait
#' association and estimation of genetic correlation using the matrix-variate
#' linear mixed model. Genetics, May 2015, Vol.200-1, p. 59-68.
#' @references Korte et al. (2012). A mixed-model approach for genome-wide
#' association studies of correlated traits in structured populations.
#' Nature Genetics, 44(9), 1066–1071. https://doi.org/10.1038/ng.2376
#' @references Millet et al. (2016). Genome-wide analysis of yield in Europe:
#' allelic effects as functions of drought and heat scenarios. Plant Physiology,
#' pp.00621.2016. https://doi.org/10.1104/pp.16.00621
#' @references Thoen et al. (2016). Genetic architecture of plant stress
#' resistance: multi-trait genome-wide association mapping. New Phytologist,
#' 213(3), 1346–1362. https://doi.org/10.1111/nph.14220
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409.
#'
#' @importFrom data.table :=
#' @importFrom methods as
#'
#' @export
runMultiTraitGwas <- function(gData,
                              environments = NULL,
                              covar = NULL,
                              snpCov = NULL,
                              kin = NULL,
                              kinshipMethod = c("astle", "GRM", "IBS",
                                                "vanRaden"),
                              GLSMethod = c("single", "multi"),
                              subsetMarkers = FALSE,
                              markerSubset = "",
                              MAF = 0.01,
                              fitVarComp = TRUE,
                              covModel = c("unst", "pw", "fa"),
                              VeDiag = TRUE,
                              tolerance = 1e-6,
                              maxIter = 2e5,
                              maxDiag = 1e4,
                              mG = 1,
                              mE = 1,
                              CmHet = TRUE,
                              DmHet = TRUE,
                              stopIfDecreasing = TRUE,
                              Vg = NULL,
                              Ve = NULL,
                              reduceK = FALSE,
                              nPca = NULL,
                              estCom = FALSE,
                              parallel = FALSE,
                              nCores = NULL) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers)
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
  chkCovar(covar, gData)
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) {
    covar <- colnames(gData$covar)[covar]
  }
  chkSnpCov(snpCov, gData)
  ## SNPs with MAF == 0 always have to be removed to prevent creation of
  ## singular matrices.
  chkNum(MAF, min = 0, max = 1)
  MAF <- max(MAF, 1e-6)
  ## If environments is null set environments to only environment in pheno.
  if (is.null(environments)) {
    environments <- 1
  }
  GLSMethod <- match.arg(GLSMethod)
  covModel <- match.arg(covModel)
  if (covModel == "fa") {
    chkNum(tolerance, min = 0)
    chkNum(maxIter, min = 1)
    chkNum(maxDiag, min = 0)
    chkNum(mG, min = 1)
    chkNum(mE, min = 1)
  }
  chkKin(kin, gData, GLSMethod)
  kinshipMethod <- match.arg(kinshipMethod)
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
                 "included in column names of pheno.\n"))
    }
    if (is.null(colnames(Ve)) || is.null(rownames(Ve)) ||
        any(colnames(Ve) != rownames(Ve)) ||
        !all(colnames(Ve) %in% colnames(gData$pheno[[1]])[-1])) {
      stop(paste("Column names and rownames of Ve should be identical and",
                 "included in column names of pheno.\n"))
    }
    Vg <- Vg[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    Ve <- Ve[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    colnames(Vg) <- rownames(Vg) <- colnames(Ve) <- rownames(Ve) <- NULL
  }
  if (reduceK) {
    chkNum(nPca, min = 1)
  }
  markers <- gData$markers
  map <- gData$map
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
  #X <- cbind(rep(1, nrow(phEnv)), as(as.matrix(phEnv[covEnv]), "dgeMatrix"))
  X <- cbind(rep(1, nrow(phEnv)), as.matrix(phEnv[covEnv]))
  rownames(X) <- phEnv$genotype
  ## Add snpCovariates to X
  if (!is.null(snpCov)) {
    if (ncol(X) == length(snpCov)) {
      XRed <- matrix(nrow = nrow(X), ncol = 0, dimnames = list(rownames(X)))
    } else {
      XRed <- X[, 1:(ncol(X) - length(snpCov)), drop = FALSE]
    }
  }
  ## Construct Y from pheno data in gData.
  Y <- phEnv[, !colnames(phEnv) %in% covEnv]
  rownames(Y) <- Y[["genotype"]]
  Y <- as.matrix(Y[, -which(colnames(Y) == "genotype")])
  if (anyNA(Y)) {
    stop("Phenotypic data cannot contain any missing values.\n")
  }
  if (GLSMethod == "single") {
    ## Compute kinship matrix.
    K <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
                    markers = markersRed, kinshipMethod = kinshipMethod)
    K <- K[rownames(K) %in% rownames(Y), colnames(K) %in% rownames(Y)]
    if (reduceK) {
      K <- reduceKinship(K = K, nPca = nPca)
    }
    Y <- Y[rownames(Y) %in% rownames(K), ]
    X <- X[rownames(X) %in% rownames(K), , drop = FALSE]
  } else if (GLSMethod == "multi") {
    ## Compute kinship matrices per chromosome. Only needs to be done once.
    chrs <- unique(mapRed$chr[rownames(mapRed) %in% colnames(markersRed)])
    KChr <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
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
  ## fit variance components.
  if (fitVarComp) {
    if (GLSMethod == "single") {
      if (covModel == "unst") {
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
      } else if (covModel == "pw") {
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
      } else if (covModel == "fa") {
        ## FA models.
        ## Including snpCovariates.
        varComp <- EMFA(y = Y, k = as.matrix(K),
                        size_param_x = X, maxIter = maxIter,
                        tolerance = tolerance, mG = mG, mE = mE, cmHet = CmHet,
                        dmHet = DmHet, maxDiag = maxDiag,
                        stopIfDecreasing = stopIfDecreasing)
        if (!is.null(snpCov)) {
          ## Without snpCovariates.
          varCompRed <- EMFA(y = Y, k = as.matrix(K),
                             size_param_x = XRed, maxIter = maxIter,
                             tolerance = tolerance, mG = mG, mE = mE,
                             cmHet = TRUE, dmHet = TRUE, maxDiag = maxDiag,
                             stopIfDecreasing = stopIfDecreasing)
        }
      }
      Vg <- varComp$Vg
      Ve <- varComp$Ve
      rownames(Vg) <- colnames(Vg) <- rownames(Ve) <- colnames(Ve) <-
        colnames(Y)
      if (!is.null(snpCov)) {
        VgRed <- varCompRed$Vg
        VeRed <- varCompRed$Ve
        rownames(VgRed) <- colnames(VgRed) <- rownames(VeRed) <-
          colnames(VeRed) <- colnames(Y)
      }
    } else if (GLSMethod == "multi") {
      if (covModel == "unst") {
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
      } else if (covModel == "pw") {
        ## Unstructured (pairwise) models.
        ## Sommer always adds an intercept so remove it from X.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          covPW(Y = Y, K = KChr[[which(chrs == chr)]],
                X = if (ncol(X) == 1) {
                  NULL
                } else {
                  X[, -1, drop = FALSE]
                }, fixDiag = FALSE, corMat = FALSE, parallel = parallel)
        }, simplify = FALSE)
        if (!is.null(snpCov)) {
          ## Sommer always adds an intercept so remove it from XRed.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            covPW(Y = Y, K = KChr[[which(chrs == chr)]],
                  X = if (ncol(XRed) == 1) {
                    NULL
                  } else {
                    XRed[, -1, drop = FALSE]
                  }, fixDiag = FALSE, corMat = FALSE, parallel = parallel)
          }, simplify = FALSE)
        }
      } else if (covModel == "fa") {
        ## FA models.
        ## Including snpCovariates.
        varComp <- sapply(X = chrs, FUN = function(chr) {
          EMFA(y = as.matrix(Y), k = as.matrix(KChr[[which(chrs == chr)]]),
               size_param_x = X, maxIter = maxIter, tolerance = tolerance,
               mG = mG, mE = mE, cmHet = CmHet, dmHet = DmHet,
               maxDiag = maxDiag, stopIfDecreasing = stopIfDecreasing)
        }, simplify = FALSE)
        if (!is.null(snpCov)) {
          ## Without snpCovariates.
          varCompRed <- sapply(X = chrs, FUN = function(chr) {
            EMFA(y = as.matrix(Y), k = as.matrix(KChr[[which(chrs == chr)]]),
                 size_param_x = XRed, maxIter = maxIter, tolerance = tolerance,
                 mG = mG, mE = mE, cmHet = CmHet, dmHet = DmHet,
                 maxDiag = maxDiag, stopIfDecreasing = stopIfDecreasing)
          }, simplify = FALSE)
        }
      }
      Vg <- setNames(lapply(X = varComp, FUN = function(Vc) {
        Vg0 <- Vc[[1]]
        rownames(Vg0) <- colnames(Vg0) <- colnames(X)
        return(Vg0)
      }), paste("chr", chrs))
      Ve <- setNames(lapply(X = varComp, FUN = function(Vc) {
        Ve0 <- Vc[[2]]
        rownames(Ve0) <- colnames(Ve0) <- colnames(X)
        return(Ve0)
      }), paste("chr", chrs))
      if (!is.null(snpCov)) {
        VgRed <- setNames(lapply(X = varCompRed, FUN = function(Vc) {
          Vg0 <- Vc[[1]]
          rownames(Vg0) <- colnames(Vg0) <- colnames(XRed)
          return(Vg0)
        }), paste("chr", chrs))
        VeRed <- setNames(lapply(X = varCompRed, FUN = function(Vc) {
          Ve0 <- Vc[[2]]
          rownames(Ve0) <- colnames(Ve0) <- colnames(XRed)
          return(Ve0)
        }), paste("chr", chrs))
      }
    } #end GLSMethod multi
  } #end varComp
  allFreq <- Matrix::colMeans(markersRed[rownames(Y),
                                         rownames(mapRed)]) / max(markersRed)
  markersRed <- markersRed[rownames(Y), ]
  ## Run GWAS.
  if (GLSMethod == "single") {
    estEffRes <- estEffTot(markers = markersRed, X = X, Y = Y, K = K,
                           XRed = XRed, Vg = Vg, Ve = Ve, snpCov = snpCov,
                           allFreq = allFreq, MAF = MAF, estCom = estCom,
                           nCores = nCores)
    list2env(estEffRes, envir = environment())
  } else if (GLSMethod == "multi") {
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
                             allFreq = allFreqChr, MAF = MAF, estCom = estCom,
                             nCores = nCores)
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
  effs <- data.table::as.data.table(effs, keep.rownames = "trait")
  effsSe <- data.table::as.data.table(effsSe, keep.rownames = "trait")
  effsLong <- data.table::melt(effs, variable.factor = FALSE,
                               variable.name = "snp", value.name = "effect",
                               id.vars = "trait")
  effsSeLong <- data.table::melt(effsSe, variable.factor = FALSE,
                                 variable.name = "snp", value.name = "effectSe",
                                 id.vars = "trait")
  effsTot <- merge(effsLong, effsSeLong)
  ## Set up a data.frame for storing results containing map info and
  ## allele frequencies.
  GWAResult <- data.table::data.table(snp = rownames(mapRed), mapRed,
                                      allFreq = allFreq, key = "snp")
  ## Merge the effects, effectsSe and pValues to the results.
  GWAResult <- merge(GWAResult, effsTot, by = "snp")
  GWAResult <- merge(GWAResult,
                     data.table::data.table(snp = names(pValues),
                                            pValue = pValues, key = "snp"))
  GWAResult[, "LOD" := -log10(GWAResult$pValue)]
  if (estCom) {
    ## Bind common effects, SE, and pvalues together.
    comDat <- data.table::data.table(snp = names(pValCom), pValCom, effsCom,
                                     effsComSe, pValQtlE)
    GWAResult <- merge(GWAResult, comDat)
  }
  ## Select and compute relevant columns.
  relCols <- c("snp", "trait", "chr", "pos", "pValue", "LOD", "effect",
               "effectSe", "allFreq",
               if (estCom) {c("pValCom", "effsCom", "effsComSe",
                              "pValQtlE")})
  data.table::setcolorder(x = GWAResult, neworder = relCols)
  data.table::setkeyv(x = GWAResult, cols = c("trait", "chr", "pos"))
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   MAF = MAF,
                   GLSMethod = GLSMethod,
                   covModel = covModel,
                   varComp = list(Vg = Vg, Ve = Ve))
  return(createGWAS(GWAResult = setNames(list(GWAResult),
                                         names(gData$pheno)[env]),
                    signSnp = NULL,
                    kin = if (GLSMethod == "single") {
                      K
                    } else {
                      KChr
                    },
                    thr = NULL,
                    GWASInfo = GWASInfo))
}
