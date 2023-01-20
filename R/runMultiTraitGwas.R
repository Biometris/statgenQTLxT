#' Perform multi-trait GWAS
#'
#' \code{runMultiTraitGwas} performs multi-trait or multi-environment Genome
#' Wide Association mapping on phenotypic and genotypic data contained in a
#' \code{gData} object.
#'
#' @section Details: runMultiTraitGwas estimates the effect of a SNP in
#' different trials or on different traits, one SNP at a time. Genetic
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
#' where \eqn{Y} is a \eqn{np \times 1} vector of phenotypic values for \eqn{n}
#' genotypes and \eqn{p} traits or trials. \eqn{x} is the \eqn{n \times 1}
#' vector of scores for the marker under consideration, and \eqn{X} the
#' \eqn{n \times q} design matrix for the other covariates. By default only a
#' trait (environment) specific intercept is included. The vector of genetic
#' background effects
#' (\eqn{\left(\begin{array}{c}G_1 \\ \vdots \\ G_p\end{array}\right)}) is
#' Gaussian with zero mean and covariance \eqn{V_g \otimes K}, where \eqn{V_g}
#' is a \eqn{p \times p} matrix of genetic (co)variances, and \eqn{K} an
#' \eqn{n \times n} kinship matrix. Similarly, the residual errors
#' (\eqn{\left(\begin{array}{c}E_1 \\ \vdots \\ E_p\end{array}\right)})
#' have covariance \eqn{V_e \otimes I_n}, for a \eqn{p \times p} matrix
#' \eqn{V_e} of residual (co)variances.
#'
#' @section Hypotheses for the SNP-effects:
#' For each SNP, the null-hypothesis \eqn{\beta_1 = \dots = \beta_p = 0} is
#' tested, using the likelihood ratio test (LRT) described in Zhou and
#' Stephens (2014). If \code{estCom = TRUE}, additional tests for a common
#' effect and for QTL x E are performed, using the parameterization
#' \eqn{\beta_j = \alpha + \alpha_j (1 \leq j \leq p)}.
#' As in Korte et al (2012), we use likelihood ratio tests, but not restricted
#' to the bivariate case. For the common effect, we fit the reduced
#' model \eqn{\beta_j = \alpha}, and test if \eqn{\alpha = 0}.
#' For QTL-by-environment interaction, we test if \eqn{\alpha_1 = \dots =
#' \alpha_p = 0}.
#'
#' @section Models for the genetic and residual covariance:
#' \eqn{V_g} and \eqn{V_e} can be provided by the user
#' (\code{fitVarComp = FALSE});
#' otherwise one of the following models is used, depending on covModel.
#' If \code{covModel = "unst"}, an unstructured model is assumed, as in Zhou and
#' Stephens (2014): \eqn{V_g} and \eqn{V_e} can be any positive-definite matrix,
#' requiring a total of \eqn{p(p + 1)/2} parameters per matrix.
#' If \code{covModel = "fa"}, a factor-analytic model is fitted using an
#' EM-algorithm, as in Millet et al (2016). \eqn{V_g} and \eqn{V_e} are assumed
#' to be of the form \eqn{W W^t + D}, where \eqn{W} is a \eqn{p \times m} matrix
#' of factor loadings and \eqn{D} a diagonal matrix with trait or environment
#' specific values. \eqn{m} is the order of the model, and the parameters
#' \code{mG} and \code{mE} specify the order used for respectively \eqn{V_g}
#' and \eqn{V_e}. \code{maxIter} sets the maximum number of iterations used
#' in the EM-algorithm.
#' Finally, if \code{covModel = "pw"}, \eqn{V_g} and \eqn{V_e} are estimated
#' 'pairwise', as in Furlotte and Eskin (2015). Looping over pairs of traits
#' or trials \eqn{1 \leq j < k \leq p},
#' \eqn{V_g[j,k] = V_g[k,j]} and \eqn{V_e[j,k] = V_e[k,j]}
#' are estimated assuming a bivariate mixed model. The diagonals of
#' \eqn{V_g} and \eqn{V_e} are fitted assuming univariate mixed models. If the
#' resulting \eqn{V_g} or \eqn{V_e} is not positive-definite, they are
#' replaced by the nearest positive-definite matrix.
#' In case \code{covModel = "unst"} or \code{"pw"} it is possible to assume
#' that \eqn{V_e} is diagonal (\code{VeDiag = TRUE})
#'
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}. The latter should not contain missing
#' values. Multi-trait or multi-environment GWAS is performed for all variables
#' in \code{pheno}.
#' @param trials A vector specifying the environment on which to run GWAS.
#' This can be either a numeric index or a character name of a list item in
#' \code{pheno}.
#' @param traits A vector of traits on which to run GWAS. These can be either
#' numeric indices or character names of columns in \code{pheno}. If
#' \code{NULL}, GWAS is run on all traits.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL}, no covariates are
#' used. An intercept is included automatically (and should not be assigned as
#' covariate). SNP-covariates should be assigned using the snpCov parameter.
#' @param snpCov An optional character vector of SNP-names to be included as
#' covariates. SNP-names should match those used in \code{gData}.
#' @param kin An optional kinship matrix or list of kinship matrices. These
#' matrices can be from the \code{matrix} class as defined in the base package
#' or from the \code{dsyMatrix} class, the class of symmetric matrices in the
#' Matrix package.\cr
#' If \code{GLSMethod} = "single" then one matrix should be provided, if
#' \code{GLSMethod} = "multi", a list of chromosome specific matrices of length
#' equal to the number of chromosomes in \code{map} in \code{gData}.\cr
#' If \code{NULL} then matrix \code{kinship} in \code{gData} is used. \cr
#' If both \code{kin} is provided and \code{gData} contains a matrix
#' \code{kinship} then \code{kin} is used.
#' @param kinshipMethod An optional character indicating the method used for
#' calculating the kinship matrix(ces). Currently "astle" (Astle and Balding,
#' 2009), "IBS", "vanRaden" (VanRaden, 2008), and "identity" are supported.
#' If a kinship matrix is supplied either in \code{gData} or in parameter
#' \code{kin}, \code{kinshipMethod} is ignored.
#' @param GLSMethod A character string indicating the method used to estimate
#' the marker effects. Either \code{single} for using a single kinship matrix,
#' or \code{multi} for using chromosome specific kinship matrices.
#' @param estCom Should the common SNP-effect model be fitted? If \code{TRUE}
#' not only the SNP-effects but also the common SNP-effect and QTL x E effect
#' are estimated.
#' @param useMAF Should the minor allele frequency be used for selecting SNPs
#' for the analysis. If \code{FALSE}, the minor allele count is used instead.
#' @param MAF The minor allele frequency (MAF) threshold used in GWAS. A
#' numerical value between 0 and 1. SNPs with MAF below this value are not taken
#' into account in the analysis, i.e. p-values and effect sizes are put to
#' missing (\code{NA}). Ignored if \code{useMAF} is \code{FALSE}.
#' @param MAC A numerical value. SNPs with minor allele count below this value
#' are not taken into account for the analysis, i.e. p-values and effect sizes
#' are set to missing (\code{NA}). Ignored if \code{useMAF} is \code{TRUE}.
#' @param fitVarComp Should the variance components be fitted? If \code{FALSE},
#' they should be supplied in \code{Vg} and \code{Ve}.
#' @param covModel A character string indicating the covariance model for the
#' genetic background (Vg) and residual effects (Ve); see details.
#' Either \code{unst} for unstructured for both Vg and
#' Ve (as in Zhou and Stephens (2014)), \code{pw} for unstructered for both Vg
#' and Ve (pairwise, as in Furlotte and Eskin (2013)) or \code{fa} for
#' factor-analytic for both Vg and Ve.\cr
#' Ignored if \code{fitVarComp} = \code{FALSE}
#' @param VeDiag Should there be environmental correlations if covModel = "unst"
#' or "pw"? If traits are measured on the same individuals, put \code{TRUE}.
#' @param maxIter An integer for the maximum number of iterations. Only used
#' when \code{covModel = "fa"}.
#' @param mG An integer. The order of the genetic part of the factor analytic
#' model. Only used when \code{covModel = "fa"}.
#' @param mE An integer. The order of the environmental part of the factor
#' analytic model. Only used when \code{covModel = "fa"}.
#' @param Vg An optional matrix with genotypic variance components. \code{Vg}
#' should have row and column names corresponding to the column names of
#' \code{gData$pheno}. It may contain additional rows and columns which will be
#' ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param Ve An optional matrix with environmental variance components.
#' \code{Ve} should have row names column names corresponding to the column
#' names of \code{gData$pheno}. It may contain additional rows and columns
#' which will be ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param genomicControl Should genomic control correction as in Devlin and
#' Roeder (1999) be applied?
#' @param thrType A character string indicating the type of threshold used for
#' the selection of candidate loci. Either \code{bonf} for using the
#' Bonferroni threshold, a LOD-threshold of \eqn{-log10(alpha/p)}, where p is
#' the number of markers and alpha can be specified in \code{alpha},
#' \code{fixed} for a self-chosen fixed LOD-threshold, specified in
#' \code{LODThr} or \code{small}, the LOD-threshold is chosen such as the SNPs
#' with the \code{nSnpLOD} smallest p-values are selected. \code{nSnpLOD} can
#' be specified.
#' @param alpha A numerical value used for calculating the LOD-threshold for
#' \code{thrType} = "bonf" and the significant p-Values for \code{thrType} =
#' "fdr".
#' @param LODThr A numerical value used as a LOD-threshold when
#' \code{thrType} = "fixed".
#' @param nSnpLOD A numerical value indicating the number of SNPs with the
#' smallest p-values that are selected when \code{thrType} = "small".
#' @param rho A numerical value used a the minimum value for SNPs to be
#' considered correlated when using \code{thrType} = "fdr".
#' @param pThr A numerical value just as the cut off value for p-Values for
#' \code{thrType} = "fdr".
#' @param sizeInclRegion An integer. Should the results for SNPs close to
#' significant SNPs be included? If so, the size of the region in centimorgan
#' or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included
#' in the region close to significant SNPs to only those SNPs that are in
#' sufficient Linkage Disequilibrium (LD) with the significant snp, where LD
#' is measured in terms of \eqn{R^2}. If for example \code{sizeInclRegion} =
#' 200000 and \code{minR2} = 0.5, then for every significant SNP also those SNPs
#' whose LD (\eqn{R^2}) with the significant SNP is at least 0.5 AND which are
#' at most 200000 away from this significant snp are included. Ignored if
#' \code{sizeInclRegion} = 0.
#' @param parallel Should the computation of variance components be done in
#' parallel? Only used if \code{covModel = "pw"}. A parallel computing
#' environment has to be setup by the user.
#' @param nCores A numerical value indicating the number of cores to be used by
#' the parallel part of the algorithm. If \code{NULL} the number of cores used
#' will be equal to the number of cores available on the machine - 1.
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
#' Nature Genetics, 44(9), 1066–1071. \doi{10.1038/ng.2376}
#' @references Millet et al. (2016). Genome-wide analysis of yield in Europe:
#' allelic effects as functions of drought and heat scenarios. Plant Physiology,
#' pp.00621.2016. \doi{10.1104/pp.16.00621}
#' @references Thoen et al. (2016). Genetic architecture of plant stress
#' resistance: multi-trait genome-wide association mapping. New Phytologist,
#' 213(3), 1346–1362. \doi{10.1111/nph.14220}
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409.
#'
#' @examples
#' ## First create a gData object.
#' ## See the vignette for a detailed description.
#' ## Here we use the gData object included in the package
#'
#' ## Run multi-trait GWAS
#' ## Use a factor analytic model to estimate variance components.
#' mtg0 <- runMultiTraitGwas(gDataDropsRestr,
#'                          trial = "Mur13W",
#'                          covModel = "fa")
#'
#' ## Plot the results.
#' ## For details on the different plots see plot.GWAS
#' plot(mtg0, plotType = "qq")
#' plot(mtg0, plotType = "manhattan")
#' plot(mtg0, plotType = "qtl", yThr = 3.5)
#'
#' ## Run multi-trait GWAS
#' ## Use a pairwise model to estimate variance components.
#' ## Estimate common effects and set a fixed threshold for significant SNPs
#' mtg1 <- runMultiTraitGwas(gDataDropsRestr,
#'                          trial = "Mur13W",
#'                          covModel = "pw",
#'                          estCom = TRUE,
#'                          thrType = "fixed",
#'                          LODThr = 3)
#'
#' ## Run multi-trait GWAS
#' ## Use an unstructured model to estimate variance components.
#' ## Identify the 5 SNPs with smallest p-values as significant SNPs.
#' ## Compute the kinship matrix using the vanRaden method.
#' \dontrun{
#' mtg2 <- runMultiTraitGwas(gDataDropsRestr,
#'                          trial = "Mur13W",
#'                          kinshipMethod = "vanRaden",
#'                          covModel = "unst",
#'                          thrType = "small",
#'                          nSnpLOD = 5)
#' }
#'
#' @importFrom data.table :=
#'
#' @export
runMultiTraitGwas <- function(gData,
                              trials = NULL,
                              traits = NULL,
                              covar = NULL,
                              snpCov = NULL,
                              kin = NULL,
                              kinshipMethod = c("astle", "IBS", "vanRaden",
                                                "identity"),
                              GLSMethod = c("single", "multi"),
                              estCom = FALSE,
                              useMAF = TRUE,
                              MAF = 0.01,
                              MAC = 10,
                              genomicControl = FALSE,
                              fitVarComp = TRUE,
                              covModel = c("unst", "pw", "fa"),
                              VeDiag = TRUE,
                              maxIter = 2e5,
                              mG = 1,
                              mE = 1,
                              Vg = NULL,
                              Ve = NULL,
                              thrType = c("bonf", "fixed", "small", "fdr"),
                              alpha = 0.05,
                              LODThr = 4,
                              nSnpLOD = 10,
                              pThr = 0.05,
                              rho = 0.4,
                              sizeInclRegion = 0,
                              minR2 = 0.5,
                              parallel = FALSE,
                              nCores = NULL) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers)
  chkTrials(trials, gData)
  if (is.null(trials) && length(gData$pheno) > 1) {
    stop("pheno contains multiple trials. Trial cannot be NULL.\n")
  }
  if (is.numeric(trials)) {
    trials <- names(gData$pheno)[trials]
  }
  ## Keep option open for extension to multiple trials.
  if (!is.null(trials)) {
    trial <- trials
  } else {
    trial <- names(gData$pheno)
  }
  chkTraits(traits, trials, gData, multi = TRUE)
  if (is.numeric(traits)) {
    ## If traits is given as numeric convert to character.
    traits <- colnames(gData$pheno[[trial]])[traits]
  } else if (is.null(traits)) {
    ## If no traits supplied extract them from pheno data.
    traits <- colnames(gData$pheno[[trial]])[-1]
  }
  nTraits <- length(traits)
  ## Restrict phenotypic data to selected traits.
  gData$pheno[[trial]] <- gData$pheno[[trial]][c("genotype", traits)]
  chkCovar(covar, gData)
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) {
    covar <- colnames(gData$covar)[covar]
  }
  chkSnpCov(snpCov, gData)
  ## SNPs with MAF == 0 always have to be removed to prevent creation of
  ## singular matrices.
  if (useMAF) {
    chkNum(MAF, min = 0, max = 1)
    MAF <- max(MAF, 1e-6)
  } else {
    chkNum(MAC, min = 0)
    MAC <- max(MAC, 1)
  }
  ## If trials is null set trials to only trial in pheno.
  if (is.null(trials)) {
    trials <- names(gData$pheno)
  }
  GLSMethod <- match.arg(GLSMethod)
  covModel <- match.arg(covModel)
  if (covModel == "fa") {
    chkNum(maxIter, min = 1)
    chkNum(mG, min = 1, max = nTraits - 1)
    chkNum(mE, min = 1, max = nTraits - 1)
  }
  if (fitVarComp && covModel == "unst") {
    if (nTraits > 5 && nTraits < 10) {
      warning("unstructured covariance models not recommended for 6 to 9 ",
              "traits. Consider using another covariance model instead.\n",
              call. = FALSE)
    } else if (nTraits > 9) {
      stop("unstructured covariance models not possible for 10 or more ",
           "traits. Try pairwise or factor analytic instead.\n")
    }
  }
  thrType <- match.arg(thrType)
  if (thrType == "bonf") {
    chkNum(alpha, min = 0)
  } else if (thrType == "fixed") {
    chkNum(LODThr, min = 0)
  } else if (thrType == "small") {
    chkNum(nSnpLOD, min = 0)
  } else if (thrType == "fdr") {
    chkNum(alpha, min = 0)
    chkNum(rho, min = 0, max = 1)
    chkNum(pThr, min = 0, max = 1)
  }
  chkNum(sizeInclRegion, min = 0)
  if (sizeInclRegion > 0) {
    chkNum(minR2, min = 0, max = 1)
  }
  chkKin(kin, gData, GLSMethod)
  kinshipMethod <- match.arg(kinshipMethod)
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
        !all(colnames(gData$pheno[[1]])[-1] %in% colnames(Vg))) {
      stop("Column names and rownames of Vg should be identical and ",
           "included in column names of pheno.\n")
    }
    if (is.null(colnames(Ve)) || is.null(rownames(Ve)) ||
        any(colnames(Ve) != rownames(Ve)) ||
        !all(colnames(gData$pheno[[1]])[-1] %in% colnames(Ve))) {
      stop("Column names and rownames of Ve should be identical and ",
           "included in column names of pheno.\n")
    }
    Vg <- Vg[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    Ve <- Ve[colnames(gData$pheno[[1]])[-1], colnames(gData$pheno[[1]])[-1]]
    colnames(Vg) <- rownames(Vg) <- colnames(Ve) <- rownames(Ve) <- NULL
    VgRed <- Vg
    VeRed <- Ve
  }
  markers <- gData$markers
  map <- gData$map
  ## Restrict map and markers to markers present in both.
  markersRed <- markers[, colnames(markers) %in% rownames(map)]
  mapRed <- map[rownames(map) %in% colnames(markers), ]
  ## Add covariates to phenotypic data.
  phExp <- expandPheno(gData = gData, trial = trial, covar = covar,
                       snpCov = snpCov)
  phTr <- phExp$phTr
  ## Restrict phenotypic data to genotypes that are actually in markers.
  phTr <- phTr[phTr[["genotype"]] %in% rownames(gData$markers), ]
  covTr <- phExp$covTr
  ## Convert pheno and covariates to format suitable for fitting var components.
  X <- cbind(rep(1, nrow(phTr)), as.matrix(phTr[covTr]))
  rownames(X) <- phTr$genotype
  ## Add snpCovariates to X.
  if (!is.null(snpCov)) {
    if (ncol(X) == length(snpCov)) {
      XRed <- matrix(nrow = nrow(X), ncol = 0, dimnames = list(rownames(X)))
    } else {
      XRed <- X[, 1:(ncol(X) - length(snpCov)), drop = FALSE]
    }
  } else {
    XRed <- NULL
  }
  ## Construct Y from pheno data in gData.
  Y <- phTr[, !colnames(phTr) %in% covTr]
  rownames(Y) <- Y[["genotype"]]
  Y <- as.matrix(Y[, -which(colnames(Y) == "genotype")])
  if (anyNA(Y)) {
    stop("Phenotypic data cannot contain any missing values.\n")
  }
  ## Compute kinship matrix (GSLMethod single)
  ## or kinship matrices per chromosome (GLSMethod multi).
  gDataRest <- gData
  if (!all(grepl(pattern = "EXT", x = colnames(gData$markers)))) {
    markersRest <- gData$markers[, !grepl(pattern = "EXT",
                                          x = colnames(gData$markers))]
    mapRest <- gData$map[!grepl(pattern = "EXT", rownames(gData$map)), ]
    gDataRest$markers <- markersRest
    gDataRest$map <- mapRest
  }
  K <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gDataRest,
                  markers = gDataRest$markers, map = gDataRest$map,
                  kinshipMethod = kinshipMethod)
  if (GLSMethod == "single") {
    Y <- Y[rownames(Y) %in% rownames(K), ]
    K <- K[rownames(K) %in% rownames(Y), colnames(K) %in% rownames(Y)]
    K <- K[rownames(Y), rownames(Y)]
    X <- X[rownames(X) %in% rownames(K), , drop = FALSE]
  } else if (GLSMethod == "multi") {
    chrs <- unique(mapRed$chr[rownames(mapRed) %in% colnames(markersRed)])
    K <- lapply(X = K, FUN = function(k) {
      k <- k[rownames(k) %in% rownames(Y), colnames(k) %in% rownames(Y)]
      k[rownames(Y), rownames(Y)]
    })
    Y <- Y[rownames(Y) %in% rownames(K[[1]]), ]
    X <- X[rownames(X) %in% rownames(K[[1]]), , drop = FALSE]
  }
  Y <- scale(Y)
  ## fit variance components.
  if (fitVarComp) {
    estVarCompRes <- estVarComp(GLSMethod = GLSMethod, covModel = covModel,
                                Y = Y, K = K, X = X, VeDiag = VeDiag,
                                snpCov = snpCov, XRed = XRed,
                                parallel = parallel, maxIter = maxIter, mG = mG,
                                mE = mE, chrs = chrs)
    list2env(estVarCompRes, envir = environment())
  }
  YScaledScale <- attr(Y, "scaled:scale")[traits]
  Y <- Y[, traits]
  maxScore <- max(markersRed)
  allFreq <- colMeans(markersRed[rownames(Y), rownames(mapRed)]) / maxScore
  markersRed <- markersRed[rownames(Y), ]
  if (!useMAF) {
    ## MAC used. Compute MAF from MAC.
    MAF <- MAC / (maxScore * nrow(Y)) - 1e-5
  }
  ## Run GWAS.
  if (GLSMethod == "single") {
    estEffRes <- estEffTot(markers = markersRed,
                           X = X, Y = Y, K = K, XRed = XRed, Vg = Vg, Ve = Ve,
                           VgRed = VgRed, VeRed = VeRed, snpCov = snpCov,
                           allFreq = allFreq, MAF = MAF, estCom = estCom,
                           nCores = nCores)
    list2env(estEffRes, envir = environment())
  } else if (GLSMethod == "multi") {
    pValues <- pValCom <- pValQtlE <- numeric()
    ## Create an empty matrix with traits as header.
    effs <- effsSe <- effsCom <- effsComSe <- t(gData$pheno[[trial]][FALSE, -1])
    for (chr in chrs) {
      mapRedChr <- mapRed[mapRed$chr == chr, ]
      markersRedChr <-
        markersRed[, colnames(markersRed) %in% rownames(mapRedChr),
                   drop = FALSE]
      allFreqChr <- colMeans(markersRedChr) / max(markersRedChr)
      snpCovChr <- snpCov[snpCov %in% colnames(markersRedChr)]
      chrNum <- which(chrs == chr)
      estEffRes <- estEffTot(markers = markersRedChr,
                             X = X, Y = Y, K = K[[chrNum]], XRed = XRed,
                             Vg = Vg[[chrNum]], Ve = Ve[[chrNum]],
                             VgRed = VgRed[[chrNum]], VeRed = VeRed[[chrNum]],
                             snpCov = snpCovChr, allFreq = allFreqChr,
                             MAF = MAF, estCom = estCom, nCores = nCores)
      pValues <- c(pValues, estEffRes$pValues)
      effs <- cbind(effs, estEffRes$effs)
      effsSe <- cbind(effsSe, estEffRes$effsSe)
      pValCom <- c(pValCom, estEffRes$pValCom)
      effsCom <- c(effsCom, estEffRes$effsCom)
      effsComSe <- c(effsComSe, estEffRes$effsComSe)
      pValQtlE <- c(pValQtlE, estEffRes$pValQtlE)
    }
  }
  ## Rescale effects and standard errors.
  effs <- effs * YScaledScale
  effsSe <- effsSe * YScaledScale
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
  ## Set up a data.table for storing results containing map info and
  ## allele frequencies.
  GWAResult <- data.table::data.table(snp = rownames(mapRed), mapRed,
                                      allFreq = allFreq, key = "snp")
  ## Merge the effects, effectsSe and pValues to the results.
  GWAResult <- merge(GWAResult, effsTot, by = "snp")
  GWAResult <- merge(GWAResult,
                     data.table::data.table(snp = names(pValues),
                                            pValue = pValues, key = "snp"))
  if (estCom) {
    ## Bind common effects, SE, and pvalues together.
    comDat <- data.table::data.table(snp = names(pValCom), pValCom, effsCom,
                                     effsComSe, pValQtlE)
    GWAResult <- merge(GWAResult, comDat)
  }
  ## Calculate the genomic inflation factor.
  GC <- statgenGWAS:::genCtrlPVals(pVals = GWAResult[GWAResult[["trait"]] == traits[1]][["pValue"]],
                                   nObs = nrow(GWAResult) / length(traits),
                                   nCov = length(covTr))
  inflationFactor <- GC$inflation
  ## Rescale p-values.
  if (genomicControl) {
    GWAResult[, "pValue" := rep(x = GC$pValues, each = length(traits))]
  }
  GWAResult[, "LOD" := -log10(GWAResult[["pValue"]])]
  ## Select relevant columns.
  relCols <- c("snp", "trait", "chr", "pos", "pValue", "LOD", "effect",
               "effectSe", "allFreq",
               if (estCom) {c("pValCom", "effsCom", "effsComSe",
                              "pValQtlE")})
  ## Reorder columns.
  data.table::setcolorder(x = GWAResult, neworder = relCols)
  ## When thrType is bonferroni or small, determine the LOD threshold.
  if (thrType == "bonf") {
    ## Compute LOD threshold using Bonferroni correction.
    LODThr <- -log10(alpha / (sum(!is.na(GWAResult[["pValue"]])) / nTraits))
  } else if (thrType == "small") {
    ## Compute LOD threshold by computing the 10log of the nSnpLOD item
    ## of ordered p values.
    LODThr <- sort(na.omit(GWAResult[GWAResult[["trait"]] == traits[1], ][["LOD"]]),
                   decreasing = TRUE)[nSnpLOD]
  } else if (thrType == "fdr") {
    LODThr <- NA_real_
  }
  ## Select the SNPs whose LOD-scores are above the threshold.
  maxScore <- max(markersRed)
  traits <- colnames(Y)
  LODThrTr <- setNames(rep(LODThr, length(traits)), traits)
  if (thrType == "fdr") {
    signSnpTr <- lapply(X = traits, FUN = function(tr) {
      extrSignSnpsFDR(GWAResult = GWAResult[GWAResult[["trait"]] == tr],
                      markers = markersRed,
                      maxScore = maxScore, pheno = phTr[, tr, drop = FALSE],
                      trait = tr, rho = rho, pThr = pThr,
                      alpha = alpha)
    })
  } else {
    signSnpTr <- lapply(X = traits, FUN = function(tr) {
      extrSignSnps(GWAResult = GWAResult[GWAResult[["trait"]] == tr],
                   LODThr = LODThr, sizeInclRegion = sizeInclRegion,
                   minR2 = minR2, map = mapRed, markers = markersRed,
                   maxScore = maxScore, pheno = phTr[, tr, drop = FALSE],
                   trait = tr)
    })
  }
  signSnp <- do.call(rbind, signSnpTr)
  ## No significant SNPs should return NULL instead of data.table().
  if (nrow(signSnp) == 0) {
    signSnp <- NULL
  }
  ## Sort columns.
  data.table::setkeyv(x = GWAResult, cols = c("trait", "chr", "pos"))
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   MAF = MAF,
                   thrType = thrType,
                   GLSMethod = GLSMethod,
                   covModel = covModel,
                   varComp = list(Vg = Vg, Ve = Ve),
                   genomicControl = genomicControl,
                   inflationFactor = inflationFactor)
  return(createGWAS(GWAResult = setNames(list(GWAResult), trials),
                    signSnp = setNames(list(signSnp), trials),
                    kin = K,
                    thr = setNames(rep(list(LODThrTr),
                                       length(trials)), trials),
                    GWASInfo = GWASInfo))
}
