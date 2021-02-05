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
#' where Y is a np x 1 vector of phenotypic values for n genotypes and p traits
#' or trials. x is the n x 1 vector of scores for the marker under
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
#' For each SNP, the null-hypothesis \eqn{\beta_1 = \dots = \beta_p = 0} is
#' tested, using the likelihood ratio test (LRT) described in Zhou and
#' Stephens (2014). If \code{estCom = TRUE}, additional tests for a common effect and
#' for QTL x E are performed, using the parameterization \eqn{\beta_j = \alpha +
#' \alpha_j (1 \leq j \leq p)}. As in Korte et al (2012), we use likelihood ratio
#' tests, but not restricted to the bivariate case. For the common effect, we
#' fit the reduced model \eqn{\beta_j = \alpha}, and test if \eqn{\alpha = 0}.
#' For QTL-by-environment interaction, we test if \eqn{\alpha_1 = \dots =
#' \alpha_p = 0}.
#'
#' @section Models for the genetic and residual covariance:
#' \eqn{V_g} and \eqn{V_e} can be provided by the user
#' (\code{fitVarComp = FALSE});
#' otherwise one of the following models is used, depending on covModel.
#' If \code{covModel = "unst"}, an unstructured model is assumed, as in Zhou and
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
#' or trials \eqn{1 \leq j < k \leq p},
#' \eqn{V_g[j,k] = V_g[k,j]} and \eqn{V_e[j,k] = V_e[k,j]}
#' are estimated assuming a bivariate mixed model. The diagonals of
#' \eqn{V_g} and \eqn{V_e} are fitted assuming univariate mixed models. If the
#' resulting \eqn{V_g} or \eqn{V_e} is not positive-definite, they are
#' replaced by the nearest positive-definite matrix.
#' In case \code{covModel = "unst"} or \code{"pw"} it is possible to assume that
#' \eqn{V_e} is diagonal (\code{VeDiag = TRUE})
#'
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}. The latter should not contain missing
#' values. Multi-trait or multi-environment GWAS is performed for all variables
#' in \code{pheno}.
#' @param trials A vector specifying the environment on which to run GWAS.
#' Thise can be either a numeric index or a character name of a list item in
#' \code{pheno}.
#' @param traits A vector of traits on which to run GWAS. These can be either
#' numeric indices or character names of columns in \code{pheno}. If \code{NULL},
#' GWAS is run on all traits.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL}, no covariates are
#' used. An intercept is included automatically (and should not be assigned as
#' covariate). SNP-covariates should be assigned using the snpCov parameter.
#' @param estCom Should the common SNP-effect model be fitted? If \code{TRUE}
#' not only the SNP-effects but also the common SNP-effect and QTL x E effect
#' are estimated.
#' @param MAF The minor allele frequency (MAF) threshold used in GWAS. A
#' numerical value between 0 and 1. SNPs with MAF below this value are not taken
#' into account in the analysis, i.e. p-values and effect sizes are put to
#' missing (\code{NA}).
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
#' @param Vg An optional matrix with genotypic variance components. \code{Vg} should
#' have row and column names corresponding to the column names of
#' \code{gData$pheno}. It may contain additional rows and columns which will be
#' ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param Ve An optional matrix with environmental variance components. \code{Ve}
#' should have row names column names corresponding to the column names of
#' \code{gData$pheno}. It may contain additional rows and columns which will be
#' ignored. Ignored if fitVarComp = \code{TRUE}.
#' @param thrType A character string indicating the type of threshold used for
#' the selection of candidate loci. Either \code{bonf} for using the
#' Bonferroni threshold, a LOD-threshold of \eqn{-log10(alpha/p)}, where p is
#' the number of markers and alpha can be specified in \code{alpha},
#' \code{fixed} for a self-chosen fixed LOD-threshold, specified in \code{LODThr}
#' or \code{small}, the LOD-threshold is chosen such as the SNPs with the
#' \code{nSnpLOD} smallest p-values are selected. \code{nSnpLOD} can be
#' specified.
#' @param alpha A numerical value used for calculating the LOD-threshold for
#' \code{thrType} = "bonf".
#' @param LODThr A numerical value used as a LOD-threshold when
#' \code{thrType} = "fixed".
#' @param nSnpLOD A numerical value indicating the number of SNPs with the
#' smallest p-values that are selected when \code{thrType} = "small".
#' @param rho A numerical value ...
#' @param pThr A numerical value ...
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
#'
#' @export
multiQTL <- function(SIM,
                     gData,
                     minCofactorProximity = 50,
                     kinshipMethod = c("astle", "IBS", "vanRaden", "identity")) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers)
  kinshipMethod <- match.arg(kinshipMethod)
  peaks <- SIM$peaks
  ## Kinship needs to be computed on observed markers.
  markersRest <- gData$markers[, !grepl(pattern = "EXT",
                                        x = colnames(gData$markers))]
  mapRest <- gData$map[!grepl(pattern = "EXT", rownames(gData$map)), ]
  ## Kinship computation requires markers on a 0-2 scale.
  markersRest <- markersRest + 1
  gDataRest <- gData
  gDataRest$markers <- markersRest
  gDataRest$map <- mapRest
  ## Compute kinship matrix (GSLMethod single)
  ## or kinship matrices per chromosome (GLSMethod multi).
  K <- computeKin(GLSMethod = "multi", gData = gDataRest,
                  markers = markersRest, map = mapRest,
                  kinshipMethod = kinshipMethod)
  trials <- names(SIM$GWAResult)
  traits <- names(gData$pheno[[1]])[names(gData$pheno[[1]]) %in%
                                      SIM$GWAResult[[trials]][["trait"]]]
  nChr <- length(unique(SIM$GWAResult[[trials]][["chr"]]))
  covModel <- SIM$GWASInfo$covModel
  MAF <- SIM$GWASInfo$MAF
  thrType <- SIM$GWASInfo$thrType
  thr <- SIM$thr[[1]][1]
  Ve <- SIM$GWASInfo$varComp$Ve
  Vg <- SIM$GWASInfo$varComp$Vg
  VeLst <- rep(list(Ve), times = nChr)
  VgLst <- rep(list(Vg), times = nChr)
  VeDiag <- all(upper.tri(Ve, diag = FALSE) == 0)
  ## Now run multi-trait GWAS.
  resGWAS <- runMultiTraitGwas(
    gData = gData, trials = trials, traits = traits,
    minCofactorProximity = minCofactorProximity,
    snpCov = unique(peaks[["snp"]]), kin = K,
    GLSMethod = "multi",
    fitVarComp = FALSE, Ve = VeLst, Vg = VgLst,
    #fitVarComp = TRUE, covModel = covModel,
    VeDiag = VeDiag, thrType = thrType, LODThr = thr)
  ## Add parent info to result.
  par1 <- attr(gData, which = "parents")[1]
  par2 <- attr(gData, which = "parents")[2]
  resGWAS$GWAResult[[1]][["highValueAllele"]] <-
    ifelse(resGWAS$GWAResult[[1]][["effect"]] > 0, par1, par2)
  resGWAS$signSnp[[1]][["highValueAllele"]] <-
    ifelse(resGWAS$signSnp[[1]][["effect"]] > 0, par1, par2)
  ## Get the peaks found.
  peaks <- getPeaks(resGWAS)
  res <- resGWAS
  res$peaks <- resGWAS$signSnp[[1]][snp %in% peaks, ]
  res$signSnp[[1]][!snp %in% peaks, "snpStatus"] <-
    "within LD of significant SNP"
  class(res) <- c("multiQTL", class(res))
  attr(res, which = "parents") <- attr(gData, which = "parents")
  return(res)
}

