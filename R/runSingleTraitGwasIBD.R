#' Perform Single Trait GWAS for IBD based QTL mapping.
#'
#' \code{runSingleTraitGwas} performs a single-trait Genome Wide Association
#' Study (GWAS) on phenotypic and genotypic data contained in a \code{gData}
#' object. A covariance matrix is computed using the EMMA algorithm (Kang et
#' al., 2008) or the Newton-Raphson algorithm (Tunnicliffe, 1989) in the
#' \code{sommer} package. Then a Generalized Least Squares (GLS) method is used
#' for estimating the marker effects and corresponding p-values. This is done
#' using either one kinship matrix for all chromosomes or different chromosome
#' specific kinship matrices for each chromosome. Significant SNPs are selected
#' based on a user defined threshold.
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param ref A numerical value indicating the allele to be used as reference
#' allele.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic
#' relatedness in genetic association studies, Stat. Sci., November 2009,
#' Vol. 24, no. 4, p. 451–471.
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association
#' studies. Biometrics, December 1999, Vol. 55(4), p. 997-1004.
#' @references Kang et al. (2008) Efficient Control of Population Structure in
#' Model Organism Association Mapping. Genetics, March 2008, Vol. 178, no. 3,
#' p. 1709-1723.
#' @references Rincent et al. (2014) Recovering power in association mapping
#' panels with variable levels of linkage disequilibrium. Genetics, May 2014.
#' Vol. 197. p. 375–387.
#' @references Segura et al. (2012) An efficient multi-locus mixed-model
#' approach for genome-wide association studies in structured populations.
#' Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association
#' mapping. Heredity, February 2010, Vol. 105, p. 333–340.
#' @references Tunnicliffe W. (1989) On the use of marginal likelihood in time
#' series model estimation. JRSS, Vol.51(1), p.15-27.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. J Dairy Sci, November 2008, Vol. 91 p. 4414–4423.
#'
#' @seealso \code{\link{GWAS}}, \code{\link{kinship}},
#' \code{\link{runSingleTraitGwas}}
#'
#' @import stats
#'
#' @export
runSingleTraitGwasIBD <- function(gData,
                                  traits = NULL,
                                  environments = NULL,
                                  covar = NULL,
                                  snpCov = NULL,
                                  kin = NULL,
                                  kinshipMethod = c("multiAllKin"),
                                  remlAlgo = c("EMMA", "NR"),
                                  GLSMethod = c("single", "multi"),
                                  useMAF = TRUE,
                                  MAF = 0.01,
                                  MAC = 10,
                                  genomicControl = FALSE,
                                  thrType = c("bonf", "fixed", "small"),
                                  alpha = 0.05 ,
                                  LODThr = 4,
                                  nSnpLOD = 10,
                                  sizeInclRegion = 0,
                                  minR2 = 0.5,
                                  ref = 1) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers, dim = 3)
  chkEnvs(environments, gData)
  ## If environments is null set environments to all environments in pheno.
  if (is.null(environments)) {
    environments <- 1:length(gData$pheno)
  }
  chkTraits(traits, environments, gData)
  chkCovar(covar, gData)
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) {
    covar <- colnames(gData$covar)[covar]
  }
  chkSnpCov(snpCov, gData)
  GLSMethod <- match.arg(GLSMethod)
  chkKin(kin, gData, GLSMethod)
  kinshipMethod <- match.arg(kinshipMethod)
  remlAlgo <- match.arg(remlAlgo)
  chkNum(sizeInclRegion, min = 0)
  if (sizeInclRegion > 0) {
    chkNum(minR2, min = 0, max = 1)
  }
  if (useMAF) {
    chkNum(MAF, min = 0, max = 1)
    MAF <- max(MAF, 1e-6)
  } else {
    chkNum(MAC, min = 0)
    MAC <- max(MAC, 1)
  }
  thrType <- match.arg(thrType)
  if (thrType == "bonf") {
    chkNum(alpha, min = 0)
  } else if (thrType == "fixed") {
    chkNum(LODThr, min = 0)
  } else if (thrType == "small") {
    chkNum(nSnpLOD, min = 0)
  }
  chkNum(ref, min = 1, max = dim(gData$markers)[3])
  if (GLSMethod == "single") {
    ## Compute kinship matrix.
    K <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
                    markers = gData$markers, map = gData$map,
                    kinshipMethod = kinshipMethod)
  } else if (GLSMethod == "multi") {
    ## Compute kinship matrices per chromosome. Only needs to be done once.
    KChr <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
                       markers = gData$markers, map = gData$map,
                       kinshipMethod = kinshipMethod)
  }
  ## Compute max value in markers.
  ###maxScore <- min(max(gData$markers, na.rm = TRUE), 2)
  ## Define data.frames for total results.
  GWATot <- signSnpTot <- varCompTot <- LODThrTot <- inflationFactorTot <-
    setNames(vector(mode = "list", length = length(environments)),
             names(gData$pheno)[environments])
  for (env in environments) {
    ## Add covariates to phenotypic data.
    phExp <- expandPheno(gData = gData, env = env, covar = covar,
                         snpCov = snpCov, ref = ref)
    phEnv <- phExp$phEnv
    covEnv <- phExp$covEnv
    ## If traits is given as numeric convert to character.
    if (is.numeric(traits)) {
      traits <- colnames(gData$pheno[[env]])[traits]
    }
    ## If no traits supplied extract them from pheno data.
    if (is.null(traits)) {
      traits <- colnames(gData$pheno[[env]])[-1]
    }
    LODThrEnv <- inflationFactorEnv <-
      setNames(numeric(length = length(traits)), traits)
    GWATotEnv <- signSnpTotEnv <- varCompEnv <-
      setNames(vector(mode = "list", length = length(traits)), traits)
    ## Perform GWAS for all traits.
    for (trait in traits) {
      ## Select relevant columns only.
      phEnvTr <- phEnv[!is.na(phEnv[trait]) & phEnv$genotype %in%
                         rownames(gData$markers), c("genotype", trait, covEnv)]
      ## Select genotypes where trait is not missing.
      nonMiss <- unique(phEnvTr$genotype)
      nonMissRepId <- phEnvTr$genotype
      if (GLSMethod == "single") {
        kinshipRed <- K[nonMiss, nonMiss]
        chrs <- NULL
      } else if (GLSMethod == "multi") {
        chrs <- unique(gData$map$chr[rownames(gData$map) %in%
                                       colnames(gData$markers)])
      }
      ## Estimate variance components.
      vc <- estVarComp(GLSMethod = GLSMethod, remlAlgo = remlAlgo,
                       trait = trait, pheno = phEnvTr, covar = covEnv,
                       K = kinshipRed, chrs = chrs, KChr = KChr,
                       nonMiss = nonMiss, nonMissRepId = nonMissRepId)
      varCompEnv[[trait]] <- vc$varComp
      vcovMatrix <- vc$vcovMatrix
      ## Compute allele frequencies based on genotypes for which phenotypic
      ## data is available.
      markersRed <-
        gData$markers[nonMiss,
                      colnames(gData$markers) %in% rownames(gData$map), ]
      mapRed <- gData$map[rownames(gData$map) %in% colnames(gData$markers), ]
      ###allFreq <- Matrix::colMeans(markersRed, na.rm = TRUE) / maxScore
      if (!useMAF) {
        MAF <- MAC / length(nonMiss) - 1e-5
      }
      ## Determine segregating markers. Exclude snps used as covariates.
      ###segMarkers <- which(allFreq >= MAF & allFreq <= (1 - MAF))
      segMarkers <- seq_along(colnames(markersRed))
      ## Create data.frame for results.
      GWAResult <- data.frame(trait = trait, snp = rownames(mapRed), mapRed,
                              pValue = NA, LOD = NA, RLR2 = NA,
                              stringsAsFactors = FALSE)
      GWAResult <- cbind(GWAResult,
                         matrix(nrow = nrow(GWAResult),
                                ncol = length(dimnames(markersRed)[[3]]) - 1,
                                dimnames =
                                  list(NULL,
                                       dimnames(markersRed)[[3]][-ref])))
      ## Define single column matrix with trait non missing values.
      y <- phEnvTr[which(phEnvTr$genotype %in% nonMiss), trait]
      if (GLSMethod == "single") {
        ## Exclude snpCovariates from segregating markers.
        exclude <- exclMarkers(snpCov = snpCov, markers = markersRed,
                               allFreq = NULL, ref = ref)
        ## The following is based on the genotypes, not the replicates:
        X <- markersRed[nonMissRepId, setdiff(segMarkers, exclude), ]
        if (length(covEnv) == 0) {
          Z <- NULL
        } else {
          ## Define covariate matrix Z.
          Z <- as.matrix(phEnvTr[which(phEnvTr$genotype %in% nonMiss), covEnv])
        }
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLSIBD(y = y, X = X, Sigma = as.matrix(vcovMatrix),
                                covs = Z, ref = ref)
        GWAResult[setdiff(segMarkers, exclude),
                  c("pValue", "RLR2", dimnames(markersRed)[[3]][-ref])] <-
          GLSResult
        ## Compute p-values and effects for snpCovariates using fastGLS.
        for (snpCovariate in snpCov) {
          GLSResultSnpCov <-
            fastGLSIBD(y = y, X = markersRed[nonMissRepId, snpCovariate, ,
                                             drop = FALSE],
                       Sigma = as.matrix(vcovMatrix),
                       covs = Z[, !grepl(pattern = paste0(snpCovariate, "_"),
                                         x = colnames(Z)), drop = FALSE],
                       ref = ref)
          GWAResult[snpCovariate,
                    c("pValue", "RLR2", dimnames(markersRed)[[3]][-ref])] <-
            GLSResultSnpCov
        }
      } else if (GLSMethod == "multi") {
        ## Similar to GLSMethod single except using chromosome specific kinship
        ## matrices.
        for (chr in chrs) {
          mapRedChr <- mapRed[which(mapRed$chr == chr), ]
          markersRedChr <- markersRed[, which(colnames(markersRed) %in%
                                                rownames(mapRedChr)), ,
                                      drop = FALSE]
          ###allFreqChr <- Matrix::colMeans(markersRedChr, na.rm = TRUE) / maxScore
          ## Determine segregating markers. Exclude snps used as covariates.
          ###segMarkersChr <- which(allFreqChr >= MAF & allFreqChr <= (1 - MAF))
          segMarkersChr <- seq_along(colnames(markersRedChr))
          ## Exclude snpCovariates from segregating markers.
          exclude <- exclMarkers(snpCov = snpCov, markers = markersRedChr,
                                 allFreq = NULL, ref = ref)
          ## Remove excluded snps from segreg markers for current chromosome.
          segMarkersChr <- setdiff(intersect(segMarkersChr,
                                             which(mapRedChr$chr == chr)),
                                   exclude)
          X <- markersRedChr[nonMissRepId, segMarkersChr, , drop = FALSE]
          if (length(covEnv) == 0) {
            Z <- NULL
          } else {
            ## Define covariate matrix Z.
            Z <- as.matrix(phEnvTr[which(phEnvTr$genotype %in% nonMiss),
                                   covEnv])
          }
          GLSResult <-
            fastGLSIBD(y = y, X = X,
                       Sigma = as.matrix(vcovMatrix[[which(chrs == chr)]]),
                       covs = Z, ref = ref)
          GWAResult[colnames(markersRedChr)[segMarkersChr],
                    c("pValue", "RLR2", dimnames(markersRedChr)[[3]][-ref])] <-
            GLSResult
          ## Compute pvalues and effects for snpCovariates using fastGLS.
          for (snpCovariate in intersect(snpCov, colnames(markersRedChr))) {
            GLSResultSnpCov <-
              fastGLSIBD(y = y, X = markersRedChr[nonMissRepId, snpCovariate, ,
                                                   drop = FALSE],
                         Sigma = as.matrix(vcovMatrix[[which(chrs == chr)]]),
                         covs = Z[, !grepl(pattern = paste0(snpCovariate, "_"),
                                           x = colnames(Z)), drop = FALSE],
                         ref = ref)
            GWAResult[snpCovariate,
                      c("pValue", "RLR2", dimnames(markersRed)[[3]][-ref])] <-
              GLSResultSnpCov
          }
        }
      }
      ## Calculate the genomic inflation factor.
      GC <- genCtrlPVals(pVals = GWAResult$pValue, nObs = length(nonMiss),
                         nCov = length(covEnv))
      inflationFactorEnv[trait] <- GC$inflation
      ## Rescale p-values.
      if (genomicControl) {
        GWAResult$pValue <- GC$pValues
      }
      ## Compute LOD score.
      GWAResult$LOD <- -log10(GWAResult$pValue)
      ## Add gene information if available.
      if (!is.null(gData$genes)) {
        GWAResult <- cbind(GWAResult, gene1 = gData$genes$gene1,
                           gene2 = gData$genes$gene2)
      }
      ## When thrType is bonferroin or small, determine the LOD threshold.
      if (thrType == "bonf") {
        ## Compute LOD threshold using Bonferroni correction.
        LODThr <- -log10(alpha / sum(!is.na(GWAResult$pValue)))
      } else if (thrType == "small") {
        ## Compute LOD threshold by computing the 10log of the nSnpLOD item
        ## of ordered p values.
        LODThr <- sort(na.omit(GWAResult$LOD), decreasing = TRUE)[nSnpLOD]
      }
      LODThrEnv[trait] <- LODThr
      ## Select the SNPs whose LOD-scores is above the threshold
      signSnpTotEnv[[trait]] <-
        extrSignSnps(GWAResult = GWAResult, LODThr = LODThr,
                     sizeInclRegion = sizeInclRegion, minR2 = minR2,
                     map = mapRed, markers = markersRed,
                     maxScore = NULL, pheno = phEnvTr, trait = trait)
      GWATotEnv[[trait]] <- GWAResult
    } # end for (trait in traits)
    ## Bind data together for results and significant SNPs.
    GWATot[[match(env, environments)]] <- dfBind(GWATotEnv)
    signSnpTot[[match(env, environments)]] <- dfBind(signSnpTotEnv)
    ## No significant SNPs should return NULL instead of data.frame().
    signSnpTot <- lapply(signSnpTot, FUN = function(x) {
      if (is.null(x) || nrow(x) == 0) NULL else x
    })
    varCompTot[[match(env, environments)]] <- varCompEnv
    LODThrTot[[match(env, environments)]] <- LODThrEnv
    inflationFactorTot[[match(env, environments)]] <- inflationFactorEnv
  } # end for (environment in environments)
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   remlAlgo = remlAlgo,
                   thrType = thrType,
                   MAF = MAF,
                   GLSMethod = GLSMethod,
                   varComp = varCompTot,
                   genomicControl = genomicControl,
                   inflationFactor = inflationFactorTot)
  return(createGWAS(GWAResult = GWATot,
                    signSnp = signSnpTot,
                    kin = if (GLSMethod == "single") {
                      K
                    } else {
                      KChr
                    },
                    thr = LODThrTot,
                    GWASInfo = GWASInfo))
}
