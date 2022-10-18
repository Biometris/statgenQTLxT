# #' Perform multi-trait GWAS
# #'
# #' \code{runMultiTraitGwas} performs multi-trait or multi-environment Genome
# #' Wide Association mapping on phenotypic and genotypic data contained in a
# #' \code{gData} object.
# #'
# #' @section Details: runMultiTraitGwas estimates the effect of a SNP in
# #' different trials or on different traits, one SNP at a time. Genetic
# #' and residual covariances are fitted only once, for a model without SNPs.
# #' Following the diagonalization scheme of Zhou and Stephens (2014), the
# #' following model is fit
# #'
# #' \eqn{Y = \left(\begin{array}{c} Y_1 \\ \vdots \\ Y_p\end{array}\right) =
# #' \left(\begin{array}{c} X_1\gamma_1 \\ \vdots \\ X_p\gamma_p\end{array}\right)
# #' + \left(\begin{array}{c} x_1\beta_1 \\ \vdots \\ x_p\beta_p\end{array}\right)
# #' + \left(\begin{array}{c} G_1 \\ \vdots \\ G_p\end{array}\right) +
# #' \left(\begin{array}{c} E_1 \\ \vdots \\ E_p\end{array}\right)}
# #'
# #' where Y is a np x 1 vector of phenotypic values for n genotypes and p traits
# #' or trials. x is the n x 1 vector of scores for the marker under
# #' consideration, and X the n x q design matrix for the other covariates. By
# #' default only a trait (environment) specific intercept is included. The vector
# #' of genetic background effects
# #' (\eqn{\left(\begin{array}{c}G_1 \\ \vdots \\ G_p\end{array}\right)}) is
# #' Gaussian with
# #' zero mean and covariance \eqn{V_g \otimes K}, where \eqn{V_g} is a p x p
# #' matrix of genetic (co)variances, and K an n x n kinship matrix. Similarly,
# #' the residual errors
# #' (\eqn{\left(\begin{array}{c}E_1 \\ \vdots \\ E_p\end{array}\right)})
# #' have covariance
# #' \eqn{V_e \otimes I_n}, for a p x p matrix \eqn{V_e} of residual
# #' (co)variances.
# #'
# #' @section Hypotheses for the SNP-effects:
# #' For each SNP, the null-hypothesis \eqn{\beta_1 = \dots = \beta_p = 0} is
# #' tested, using the likelihood ratio test (LRT) described in Zhou and
# #' Stephens (2014). If \code{estCom = TRUE}, additional tests for a common effect and
# #' for QTL x E are performed, using the parameterization \eqn{\beta_j = \alpha +
# #' \alpha_j (1 \leq j \leq p)}. As in Korte et al (2012), we use likelihood ratio
# #' tests, but not restricted to the bivariate case. For the common effect, we
# #' fit the reduced model \eqn{\beta_j = \alpha}, and test if \eqn{\alpha = 0}.
# #' For QTL-by-environment interaction, we test if \eqn{\alpha_1 = \dots =
# #' \alpha_p = 0}.
# #'
# #' @section Models for the genetic and residual covariance:
# #' \eqn{V_g} and \eqn{V_e} can be provided by the user
# #' (\code{fitVarComp = FALSE});
# #' otherwise one of the following models is used, depending on covModel.
# #' If \code{covModel = "unst"}, an unstructured model is assumed, as in Zhou and
# #' Stephens (2014): \eqn{V_g} and \eqn{V_e} can be any positive-definite matrix,
# #' requiring a total of p(p+1)/2 parameters per matrix.
# #' If \code{covModel = "fa"}, a factor-analytic model is fitted using an
# #' EM-algorithm, as in Millet et al (2016). \eqn{V_g} and \eqn{V_e} are assumed
# #' to be of the form \eqn{W W^t + D}, where W is a p x m matrix of factor
# #' loadings and D a diagonal matrix with trait or environment specific values.
# #' m is the order of the model, and the parameters \code{mG} and \code{mE}
# #' specify the order used for respectively \eqn{V_g} and \eqn{V_e}.
# #' \code{maxIter} sets the maximum number of iterations used in the
# #' EM-algorithm.
# #' Finally, if \code{covModel = "pw"}, \eqn{V_g} and \eqn{V_e} are estimated
# #' 'pairwise', as in Furlotte and Eskin (2015). Looping over pairs of traits
# #' or trials \eqn{1 \leq j < k \leq p},
# #' \eqn{V_g[j,k] = V_g[k,j]} and \eqn{V_e[j,k] = V_e[k,j]}
# #' are estimated assuming a bivariate mixed model. The diagonals of
# #' \eqn{V_g} and \eqn{V_e} are fitted assuming univariate mixed models. If the
# #' resulting \eqn{V_g} or \eqn{V_e} is not positive-definite, they are
# #' replaced by the nearest positive-definite matrix.
# #' In case \code{covModel = "unst"} or \code{"pw"} it is possible to assume that
# #' \eqn{V_e} is diagonal (\code{VeDiag = TRUE})
# #'
# #' @inheritParams SIM
# #'
# #' @param SIM An object of class \code{\link{SIM}}
# #' @param gData An object of class \code{gData} containing at least \code{map},
# #' \code{markers} and \code{pheno}. The latter should not contain missing
# #' values. Multi-trait or multi-environment GWAS is performed for all variables
# #' in \code{pheno}.
# #' @param minCofactorProximity A numerical value ...
# #'
# #' @return An object of class \code{\link{SIM}}.
# #'
# #' @references Dahl et al. (2013). Network inference in matrix-variate Gaussian
# #' models with non-independent noise. arXiv preprint arXiv:1312.1622.
# #' @references Furlotte, N.A. and Eskin, E. (2015). Efficient multiple-trait
# #' association and estimation of genetic correlation using the matrix-variate
# #' linear mixed model. Genetics, May 2015, Vol.200-1, p. 59-68.
# #' @references Korte et al. (2012). A mixed-model approach for genome-wide
# #' association studies of correlated traits in structured populations.
# #' Nature Genetics, 44(9), 1066–1071. https://doi.org/10.1038/ng.2376
# #' @references Millet et al. (2016). Genome-wide analysis of yield in Europe:
# #' allelic effects as functions of drought and heat scenarios. Plant Physiology,
# #' pp.00621.2016. https://doi.org/10.1104/pp.16.00621
# #' @references Thoen et al. (2016). Genetic architecture of plant stress
# #' resistance: multi-trait genome-wide association mapping. New Phytologist,
# #' 213(3), 1346–1362. https://doi.org/10.1111/nph.14220
# #' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
# #' mixed model algorithms for genome-wide association studies. Nature Methods,
# #' February 2014, Vol. 11, p. 407–409.
# #'
# #' @importFrom data.table :=
# #'
# #' @export
# multiQTL <- function(SIM,
#                      gData,
#                      minCofactorProximity = 50,
#                      kinshipMethod = c("astle", "IBS", "vanRaden", "identity")) {
#   ## Checks.
#   chkGData(gData)
#   chkMarkers(gData$markers)
#   kinshipMethod <- match.arg(kinshipMethod)
#   peaks <- SIM$peaks
#   ## Kinship needs to be computed on observed markers.
#   markersRest <- gData$markers[, !grepl(pattern = "EXT",
#                                         x = colnames(gData$markers))]
#   mapRest <- gData$map[!grepl(pattern = "EXT", rownames(gData$map)), ]
#   ## Kinship computation requires markers on a 0-2 scale.
#   markersRest <- markersRest + 1
#   gDataRest <- gData
#   gDataRest$markers <- markersRest
#   gDataRest$map <- mapRest
#   ## Compute kinship matrix (GSLMethod single)
#   ## or kinship matrices per chromosome (GLSMethod multi).
#   K <- computeKin(GLSMethod = "multi", gData = gDataRest,
#                   markers = markersRest, map = mapRest,
#                   kinshipMethod = kinshipMethod)
#   trials <- names(SIM$GWAResult)
#   traits <- names(gData$pheno[[1]])[names(gData$pheno[[1]]) %in%
#                                       SIM$GWAResult[[trials]][["trait"]]]
#   nChr <- length(unique(SIM$GWAResult[[trials]][["chr"]]))
#   covModel <- SIM$GWASInfo$covModel
#   MAF <- SIM$GWASInfo$MAF
#   thrType <- SIM$GWASInfo$thrType
#   thr <- SIM$thr[[1]][1]
#   Ve <- SIM$GWASInfo$varComp$Ve
#   Vg <- SIM$GWASInfo$varComp$Vg
#   VeLst <- rep(list(Ve), times = nChr)
#   VgLst <- rep(list(Vg), times = nChr)
#   VeDiag <- all(upper.tri(Ve, diag = FALSE) == 0)
#   ## Now run multi-trait GWAS.
#   resGWAS <- runMultiTraitGwas(
#     gData = gData, trials = trials, traits = traits,
#     minCofactorProximity = minCofactorProximity,
#     snpCov = unique(peaks[["snp"]]), kin = K,
#     GLSMethod = "multi",
#     fitVarComp = FALSE, Ve = VeLst, Vg = VgLst,
#     VeDiag = VeDiag, thrType = thrType, LODThr = thr)
#   ## Add parent info to result.
#   par1 <- attr(gData, which = "parents")[1]
#   par2 <- attr(gData, which = "parents")[2]
#   resGWAS$GWAResult[[1]][["highValueAllele"]] <-
#     ifelse(resGWAS$GWAResult[[1]][["effect"]] > 0, par1, par2)
#   resGWAS$signSnp[[1]][["highValueAllele"]] <-
#     ifelse(resGWAS$signSnp[[1]][["effect"]] > 0, par1, par2)
#   ## Get the peaks found.
#   peaks <- getPeaks(resGWAS)
#   res <- resGWAS
#   QTLDat <- resGWAS$signSnp[[1]][resGWAS$signSnp[[1]][["snp"]] %in% peaks, ]
#   ## Fit model to compute explained variance.
#   X <- gData$markers[, colnames(gData$markers) %in% peaks, drop = FALSE]
#   modDat <- cbind(gData$pheno[[1]], X)
#   explVar <- lapply(X = traits, FUN = function(tr) {
#     ## Construct model formula.
#     ## All QTLS as random effects.
#     modFormTr <- formula(paste(tr, "~ 1 + ",
#                              paste("(1|", peaks, ")", collapse = "+")))
#     modTr <- lme4::lmer(modFormTr, data = modDat)
#     vcTr <- as.data.frame(lme4::VarCorr(modTr))
#     ## Compute explained variance.
#     vcTr[["pcts"]] <- vcTr[, "vcov"] / sum(vcTr[["vcov"]])
#     vcTr[["trait"]] <- tr
#     return(vcTr[, c("trait", "grp", "pcts")])
#   })
#   ## Bind all variances together.
#   explVar <- do.call(rbind, explVar)
#   ## Merge to QTLS.
#   QTLDat <- merge(QTLDat, explVar, by.x = c("trait", "snp"),
#                   by.y = c("trait", "grp"), sort = FALSE)
#   ## Rename columns to match other functions.
#   QTLDat[["propSnpVar"]] <- QTLDat[["pcts"]]
#   QTLDat[["pcts"]]<- NULL
#   res$peaks <- QTLDat
#   res$signSnp[[1]][!res$signSnp[[1]][["snp"]] %in% peaks, "snpStatus"] <-
#     "within LD of significant SNP"
#   trtMeans <- sapply(X = gData$pheno[[1]][, -1, drop = FALSE],
#                      FUN = mean, na.rm = TRUE)
#   class(res) <- c("multiQTL", class(res))
#   attr(res, which = "parents") <- attr(gData, which = "parents")
#   attr(res, which = "trtMeans") <- trtMeans
#   return(res)
# }

