#' Perform IBD based QTL mapping for single trait.
#'
#' \code{runSingleTraitGwasIBD} QTL mapping based on IBD probabilities.
#'
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}.
#' @param traits A vector of traits on which to run GWAS. These can be either
#' numeric indices or character names of columns in \code{pheno}. If \code{NULL},
#' GWAS is run on all traits.
#' @param trials A vector of trials on which to run GWAS. These can
#' be either numeric indices or character names of list items in \code{pheno}.
#' If \code{NULL}, GWAS is run for all trials. GWAS is run for the
#' selected trials in sequential order.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL} no covariates are
#' used.
#' @param snpCov An optional character vector of snps to be included as
#' covariates.
#' @param computeVarComp Should the variance components be computed in the
#' algorithm. If \code{FALSE} specify \code{varCovMatrix}.
#' @param varCovMatrix The variance covariance matrix used in the algorithm.
#' Providing this skips the computation of this matrix by the algorithm.
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
#' 2009), "IBS" and "vanRaden" (VanRaden, 2008) are supported. If a
#' kinship matrix is supplied either in \code{gData} or in parameter \code{kin},
#' \code{kinshipMethod} is ignored.
#' @param remlAlgo A character string indicating the algorithm used to estimate
#' the variance components. Either \code{EMMA}, for the EMMA algorithm, or
#' \code{NR}, for the Newton-Raphson algorithm.
#' @param GLSMethod A character string indicating the method used to estimate
#' the marker effects. Either \code{single} for using a single kinship matrix,
#' or \code{multi} for using chromosome specific kinship matrices.
#' @param MAF The minor allele frequency (MAF) threshold used in GWAS. A
#' numerical value between 0 and 1. SNPs with MAF below this value are not taken
#' into account in the analysis, i.e. p-values and effect sizes are put to
#' missing (\code{NA}). Ignored if \code{useMAF} is \code{FALSE}.
#' @param genomicControl Should genomic control correction as in Devlin and
#' Roeder (1999) be applied?
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
#' @param nCores A numerical value indicating the number of cores to be used by
#' the parallel part of the algorithm. If \code{NULL} the number of cores used
#' will be equal to the number of cores available on the machine - 1.
#' @param ref A numerical value indicating the founder allele to be used as
#' reference allele. Effects for the other founder allele are computed based on
#' this reference allele.
#'
#' @return An object of class \code{\link{GWAS}}.
#'
#' @seealso \code{\link{GWAS}}, \code{\link{kinship}}
#'
#' @import stats
#'
#' @export
runSingleTraitGwasIBD <- function(gData,
                                  traits = NULL,
                                  trials = NULL,
                                  covar = NULL,
                                  snpCov = NULL,
                                  computeVarComp = TRUE,
                                  kin = NULL,
                                  kinshipMethod = c("multiAllKin"),
                                  remlAlgo = c("EMMA", "NR"),
                                  GLSMethod = c("single", "multi"),
                                  varCovMatrix = NULL,
                                  MAF = 0.01,
                                  genomicControl = FALSE,
                                  thrType = c("bonf", "fixed", "small"),
                                  alpha = 0.05 ,
                                  LODThr = 4,
                                  nSnpLOD = 10,
                                  sizeInclRegion = 0,
                                  minR2 = 0.5,
                                  ref = 1,
                                  nCores = NULL) {
  ## Checks.
  chkGData(gData)
  chkMarkers(gData$markers, dim = 3)
  chkTrials(trials, gData)
  ## If trials is null set trials to all trials in pheno.
  if (is.null(trials)) {
    trials <- 1:length(gData$pheno)
  }
  chkTraits(traits, trials, gData, multi = TRUE)
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
  chkNum(MAF, min = 0, max = 1)
  MAF <- max(MAF, 1e-6)
  thrType <- match.arg(thrType)
  if (thrType == "bonf") {
    chkNum(alpha, min = 0)
  } else if (thrType == "fixed") {
    chkNum(LODThr, min = 0)
  } else if (thrType == "small") {
    chkNum(nSnpLOD, min = 0)
  }
  chkNum(ref, min = 1, max = dim(gData$markers)[3])
  ## Kinship matrices are only needed for computation of variance components.
  ## If varcomp matrix is provided by user there is no need to compute them.
  if (computeVarComp) {
    ## Compute kinship matrix.
    K <- computeKin(GLSMethod = GLSMethod, kin = kin, gData = gData,
                    markers = gData$markers, map = gData$map,
                    kinshipMethod = kinshipMethod)
  }
  ## Compute max value in markers.
  ###maxScore <- min(max(gData$markers, na.rm = TRUE), 2)
  ## Define data.frames for total results.
  GWATot <- signSnpTot <- varCompTot <- LODThrTot <- inflationFactorTot <-
    setNames(vector(mode = "list", length = length(trials)),
             names(gData$pheno)[trials])
  for (trial in trials) {
    ## Add covariates to phenotypic data.
    phExp <- expandPheno(gData = gData, trial = trial, covar = covar,
                         snpCov = snpCov)
    phTr <- phExp$phTr
    covTr <- phExp$covTr
    ## If traits is given as numeric convert to character.
    if (is.numeric(traits)) {
      traits <- colnames(gData$pheno[[trial]])[traits]
    }
    ## If no traits supplied extract them from pheno data.
    if (is.null(traits)) {
      traits <- colnames(gData$pheno[[trial]])[-1]
    }
    LODThrTr <- inflationFactorTr <-
      setNames(numeric(length = length(traits)), traits)
    GWATotTr <- signSnpTotTr <- varCompTr <-
      setNames(vector(mode = "list", length = length(traits)), traits)
    ## Perform GWAS for all traits.
    for (trait in traits) {
      ## Select relevant columns only.
      phTrTr <- phTr[!is.na(phTr[trait]) & phTr$genotype %in%
                         rownames(gData$markers), c("genotype", trait, covTr)]
      ## Select genotypes where trait is not missing.
      nonMiss <- unique(phTrTr$genotype)
      nonMissRepId <- phTrTr$genotype
      if (computeVarComp) {
        if (GLSMethod == "single") {
          kinshipRed <- K[nonMiss, nonMiss]
          chrs <- NULL
        } else if (GLSMethod == "multi") {
          chrs <- unique(gData$map$chr[rownames(gData$map) %in%
                                         colnames(gData$markers)])
        }
        ## Estimate variance components.
        vc <- estVarComp(GLSMethod = GLSMethod, remlAlgo = remlAlgo,
                         trait = trait, pheno = phTrTr, covar = covTr,
                         K = kinshipRed, chrs = chrs, nonMiss = nonMiss,
                         nonMissRepId = nonMissRepId)
        varCompTr[[trait]] <- vc$varComp
        vcovMatrix <- vc$vcovMatrix
      } else {
        vcovMatrix <- varCovMatrix
      }
      ## Compute allele frequencies based on genotypes for which phenotypic
      ## data is available.
      markersRed <-
        gData$markers[nonMiss,
                      colnames(gData$markers) %in% rownames(gData$map), ]
      mapRed <- gData$map[rownames(gData$map) %in% colnames(gData$markers), ]
      ## Determine segregating markers. Exclude snps used as covariates.
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
      y <- phTrTr[which(phTrTr$genotype %in% nonMiss), trait]
      if (GLSMethod == "single") {
        ## Exclude snpCovariates from segregating markers.
        exclude <- exclMarkers(snpCov = snpCov, markers = markersRed,
                               allFreq = NULL, ref = ref)
        ## The following is based on the genotypes, not the replicates:
        X <- markersRed[nonMissRepId, setdiff(segMarkers, exclude), ]
        if (length(covTr) == 0) {
          Z <- NULL
        } else {
          ## Define covariate matrix Z.
          Z <- as.matrix(phTrTr[which(phTrTr$genotype %in% nonMiss), covTr])
        }
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLSIBD(y = y, X = X, Sigma = as.matrix(vcovMatrix),
                                covs = Z, ref = ref, nCores = nCores)
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
                       ref = ref, nCores = nCores)
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
          ## Determine segregating markers. Exclude snps used as covariates.
          segMarkersChr <- seq_along(colnames(markersRedChr))
          ## Exclude snpCovariates from segregating markers.
          exclude <- exclMarkers(snpCov = snpCov, markers = markersRedChr,
                                 allFreq = NULL, ref = ref)
          ## Remove excluded snps from segreg markers for current chromosome.
          segMarkersChr <- setdiff(intersect(segMarkersChr,
                                             which(mapRedChr$chr == chr)),
                                   exclude)
          X <- markersRedChr[nonMissRepId, segMarkersChr, , drop = FALSE]
          if (length(covTr) == 0) {
            Z <- NULL
          } else {
            ## Define covariate matrix Z.
            Z <- as.matrix(phTrTr[which(phTrTr$genotype %in% nonMiss),
                                   covTr])
          }
          GLSResult <-
            fastGLSIBD(y = y, X = X,
                       Sigma = as.matrix(vcovMatrix[[which(chrs == chr)]]),
                       covs = Z, ref = ref, nCores = nCores)
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
                         ref = ref, nCores = nCores)
            GWAResult[snpCovariate,
                      c("pValue", "RLR2", dimnames(markersRed)[[3]][-ref])] <-
              GLSResultSnpCov
          }
        }
      }
      ## Calculate the genomic inflation factor.
      GC <- genCtrlPVals(pVals = GWAResult$pValue, nObs = length(nonMiss),
                         nCov = length(covTr))
      inflationFactorTr[trait] <- GC$inflation
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
      LODThrTr[trait] <- LODThr
      ## Select the SNPs whose LOD-scores is above the threshold
      signSnpTotTr[[trait]] <-
        extrSignSnps(GWAResult = GWAResult, LODThr = LODThr,
                     sizeInclRegion = sizeInclRegion, minR2 = minR2,
                     map = mapRed, markers = markersRed,
                     maxScore = NULL, pheno = phTrTr, trait = trait)
      GWATotTr[[trait]] <- GWAResult
    } # end for (trait in traits)
    ## Bind data together for results and significant SNPs.
    GWATot[[match(trial, trials)]] <- dfBind(GWATotTr)
    signSnpTot[[match(trial, trials)]] <- dfBind(signSnpTotTr)
    ## No significant SNPs should return NULL instead of data.frame().
    signSnpTot <- lapply(signSnpTot, FUN = function(x) {
      if (is.null(x) || nrow(x) == 0) NULL else x
    })
    varCompTot[[match(trial, trials)]] <- varCompTr
    LODThrTot[[match(trial, trials)]] <- LODThrTr
    inflationFactorTot[[match(trial, trials)]] <- inflationFactorTr
  } # end for (trial in trials)
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   remlAlgo = remlAlgo,
                   thrType = thrType,
                   MAF = MAF,
                   GLSMethod = GLSMethod,
                   varComp = varCompTot,
                   genomicControl = genomicControl,
                   inflationFactor = inflationFactorTot,
                   founders = dimnames(markersRed)[[3]][-ref])
  return(createGWAS(GWAResult = GWATot,
                    signSnp = signSnpTot,
                    kin = K,
                    thr = LODThrTot,
                    GWASInfo = GWASInfo))
}
