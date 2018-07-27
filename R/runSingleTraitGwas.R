#' Perform single-trait GWAS
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
#' @param gData An object of class \code{gData} containing at least \code{map},
#' \code{markers} and \code{pheno}.
#' @param traits A vector of traits on which to run GWAS. These can be either
#' numeric indices or character names of columns in \code{pheno}. If \code{NULL}
#' GWAS is run on all traits.
#' @param environments A vector of environments on which to run GWAS. These can
#' be either numeric indices or character names of list items in \code{pheno}.
#' If \code{NULL} GWAS is run for all environments.
#' @param covar An optional vector of covariates taken into account when
#' running GWAS. These can be either numeric indices or character names of
#' columns in \code{covar} in \code{gData}. If \code{NULL} no covariates are
#' used.
#' @param snpCov An optional character vector of snps to be included as
#' covariates.
#' @param kin An optional kinship matrix or list of kinship matrices. These
#' matrices can be from the \code{matrix} class as defined in the base package
#' or from the \code{dsyMatrix} class, the class of symmetric matrices in the
#' Matrix package.\cr
#' If \code{GLSMethod} = 1 then one matrix should be provided, if
#' \code{GLSMethod} = 2 a list of chromosome specific matrices of lenght equal
#' to the number of chromosomes in \code{map} in \code{gData}.\cr
#' If \code{NULL} then matrix \code{kinship} in \code{gData} is used. \cr
#' If both \code{kin} is provided and \code{gData} contains a matrix
#' \code{kinship} then \code{kin} is used.
#' @param kinshipMethod An optional character indicating the method used for
#' calculating the kinship matrix(ces). Currently "astle" (Astle and Balding,
#' 2009), "GRM", "IBS" and "vanRaden" (VanRaden, 2008) are supported. If a
#' kinship matrix is supplied either in \code{gData} or in parameter \code{kin}
#' \code{kinshipMethod} is ignored.
#' @param remlAlgo An integer indicating the algorithm used to estimate the
#' variance components.
#' \enumerate{
#' \item{EMMA}
#' \item{Newton-Raphson using sommer package}
#' }
#' @param GLSMethod An integer indicating the method used to estimate the marker
#' effects.
#' \enumerate{
#' \item{using one kinship matrix.}
#' \item{using chromosome specific kinship matrices.}
#' }
#' @param useMAF Should the minor allele frequency be used for selecting SNPs
#' for the analysis. If \code{FALSE} the minor allele count is used instead.
#' @param MAF A numerical value between 0 and 1. SNPs with minor allele
#' frequency below this value are not taken into account for the analysis,
#' i.e. p-values and effect sizes are put to missing (NA). Ignored if
#' \code{useMAF} is \code{FALSE}.
#' @param MAC A numerical value. SNPs with minor allele count below this value
#' are not taken into account for the analysis, i.e. p-values and effect sizes
#' are put to missing (NA). Ignored if \code{useMAF} is \code{TRUE}.
#' @param genomicControl Should genomic control correction as in Devlin and
#' Roeder (1999) be applied?
#' @param thrType An integer indicating the type of threshold used for the
#' selection of candidate loci.
#' \enumerate{
#' \item{Bonferroni; a LOD-threshold of \eqn{-log10(alpha/p)}, where p is the
#' number of markers and alpha can be specified in \code{alpha}.}
#' \item{a self-chosen LOD-threshold, specied in \code{LODThr}.}
#' \item{the LOD-threshold is chosen such the the SNPs with the \code{nSnpLOD}
#' smallest p-values are selected. \code{nSnpLOD} can be specified.}
#' }
#' @param alpha A numerical value used for calculating the LOD-threshold for
#' \code{thrType} = 1.
#' @param LODThr A numerical value used as a LOD-threshold when
#' \code{thrType} = 2.
#' @param nSnpLOD A numerical value indicating the number of SNPs with the
#' smallest p-values that are selected when \code{thrType} = 3.
#' @param sizeInclRegion An integer. Should the results for SNPs close to
#' significant SNPs be included? If so, the size of the region in centimorgan
#' or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included
#' in the region close to significant SNPs to only those SNPs that are in
#' sufficient Linkage Disequilibruim (LD) with the significant snp, where LD
#' is measured in terms of r^2. If for example \code{sizeInclRegion} = 200000
#' and \code{minR2} = 0.5, then for every significant SNP also those SNPs whose
#' LD (\eqn{r^2}) with the significant SNP is at least 0.5 AND which are at
#' most 200kb away from this significant snp are included. Ignored if
#' \code{sizeInclRegion} = 0.
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
#' @seealso \code{\link{GWAS}}, \code{\link{kinship}}
#'
#' @import stats
#'
#' @export
runSingleTraitGwas <- function(gData,
                               traits = NULL,
                               environments = NULL,
                               covar = NULL,
                               snpCov = NULL,
                               kin = NULL,
                               kinshipMethod = c("astle", "GRM", "IBS",
                                                 "vanRaden"),
                               remlAlgo = 1,
                               GLSMethod = 1,
                               useMAF = TRUE,
                               MAF = 0.01,
                               MAC = 10,
                               genomicControl = FALSE,
                               thrType = 1,
                               alpha = 0.05 ,
                               LODThr = 4,
                               nSnpLOD = 10,
                               sizeInclRegion = 0,
                               minR2 = 0.5) {
  ## Checks.
  if (missing(gData) || !is.gData(gData) || is.null(gData$map) ||
      is.null(gData$markers) || is.null(gData$pheno)) {
    stop(paste("gData should be a valid gData object with at least map,",
               "markers and pheno included.\n"))
  }
  if (!inherits(gData$markers, "Matrix")) {
    stop(paste("markers in gData should be a numerical matrix. Use",
               "recodeMarkers first for recoding.\n"))
  }
  if (anyNA(gData$markers)) {
    stop("markers contains missing values. Impute or remove these first.\n")
  }
  if (!is.null(environments) && !is.numeric(environments) &&
      !is.character(environments)) {
    stop("environments should be a numeric or character vector.\n")
  }
  if ((is.character(environments) &&
       !all(environments %in% names(gData$pheno))) ||
      (is.numeric(environments) && any(environments > length(gData$pheno)))) {
    stop("environments should be list items in pheno.\n")
  }
  ## If environments is null set environments to all environments in pheno.
  if (is.null(environments)) {
    environments <- 1:length(gData$pheno)
  }
  if (!is.null(traits) && !is.numeric(traits) && !is.character(traits)) {
    stop("traits should be a numeric or character vector.\n")
  }
  for (env in environments) {
    if ((is.character(traits) &&
         !all(traits %in% colnames(gData$pheno[[env]]))) ||
        (is.numeric(traits) &&
         (any(traits == 1) || any(traits > ncol(gData$pheno[[env]]))))) {
      stop("traits should be columns in pheno.\n")
    }
  }
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
    stop("All snpCov should be in markers.\n")
  }
  if (is.null(remlAlgo) || length(remlAlgo) > 1 || !is.numeric(remlAlgo)) {
    stop("remlAlgo should be a single numerical value.\n")
  }
  if (is.null(GLSMethod) || length(GLSMethod) > 1 || !is.numeric(GLSMethod)) {
    stop("GLSMethod should be a single numerical value.\n")
  }
  if (GLSMethod == 1 && !is.null(kin) && !(inherits(kin, "Matrix") ||
                                           is.matrix(kin))) {
    stop("kin should be a matrix.\n")
  }
  if (GLSMethod == 2 && !is.null(kin) &&
      (!is.list(kin) || !all(sapply(kin, FUN = function(k) {
        is.matrix(k) || inherits(k, "Matrix")})) ||
       length(kin) != length(unique(gData$map$chr)))) {
    stop(paste("kin should be a list of matrices of length equal to the",
               "number of chromosomes in the map.\n"))
  }
  if ((GLSMethod == 1 && is.null(gData$kinship) && is.null(kin)) ||
      GLSMethod == 2 && is.null(kin)) {
    kinshipMethod <- match.arg(kinshipMethod)
  }
  if (is.null(sizeInclRegion) || length(sizeInclRegion) > 1 ||
      !is.numeric(sizeInclRegion) || round(sizeInclRegion) != sizeInclRegion) {
    stop("sizeInclRegion should be a single integer\n")
  }
  if (sizeInclRegion > 0) {
    if (missing(minR2) || length(minR2) > 1 || !is.numeric(minR2) ||
        minR2 < 0 || minR2 > 1) {
      stop("minR2 should be a single numerical value between 0 and 1.\n")
    }
  }
  if (useMAF) {
    if (is.null(MAF) || length(MAF) > 1 || !is.numeric(MAF) ||
        MAF < 0 || MAF > 1) {
      stop("MAF should be a single numerical value between 0 and 1.\n")
    }
    if (MAF <= 1e-6) {
      MAF <- 1e-6
    }
  } else {
    if (is.null(MAC) || length(MAC) > 1 || !is.numeric(MAC)) {
      stop("MAF should be a single numerical value.\n")
    }
    if (MAC == 0) {
      MAC <- 1
    }
  }
  if (is.null(thrType) || length(thrType) > 1 || !is.numeric(thrType)) {
    stop("thrType should be a single numerical value.\n")
  }
  if (thrType == 1) {
    if (is.null(alpha) || length(alpha) > 1 || !is.numeric(alpha)) {
      stop("alpha should be a single numerical value.\n")
    } else if (thrType == 2) {
      if (is.null(LODThr) || length(LODThr) > 1 || !is.numeric(LODThr)) {
        stop("LODThr should be a single numerical value.\n")
      }
    } else if (thrType == 3) {
      if (is.null(nSnpLOD) || length(nSnpLOD) > 1 || !is.numeric(nSnpLOD)) {
        stop("nSnpLOD should be a single numerical value.\n")
      }
    }
  }
  if (!(remlAlgo %in% 1:2)) {
    remlAlgo <- 1
    warning("Invalid value for remlAlgo. remlAlgo set to 1.\n")
  }
  if (!(GLSMethod %in% 1:2)) {
    GLSMethod <- 1
    warning("Invalid value for GLSMethod. GLSMethod set to 1.\n")
  }
  if (!(thrType %in% 1:3)) {
    thrType <- 1
    warning("Invalid value for thrType. thrType set to 1.\n")
  }
  if (GLSMethod == 1) {
    ## Compute kinship matrix.
    K <- computeKin(GLSMethod = 1, kin = kin, gData = gData,
                    markers = gData$markers, kinshipMethod = kinshipMethod)
  } else if (GLSMethod == 2) {
    ## Compute kinship matrices per chromosome. Only needs to be done once.
    KChr <- computeKin(GLSMethod = 2, kin = kin, gData = gData,
                       markers = gData$markers, map = gData$map,
                       kinshipMethod = kinshipMethod)
  }
  ## Compute max value in markers
  maxScore <- min(max(gData$markers, na.rm = TRUE), 2)
  ## Define data.frames for total results.
  GWATot <- signSnpTot <- varCompTot <- LODThrTot <- inflationFactorTot <-
    setNames(vector(mode = "list", length = length(environments)),
             names(gData$pheno)[environments])
  for (env in environments) {
    ## Add covariates to phenotypic data.
    phExp <- expandPheno(gData = gData, env = env, covar = covar,
                         snpCov = snpCov)
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
      if (GLSMethod == 1) {
        kinshipRed <- K[nonMiss, nonMiss]
        chrs <- NULL
      } else if (GLSMethod == 2) {
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
      markersRed <- gData$markers[nonMiss, colnames(gData$markers) %in%
                                    rownames(gData$map)]
      mapRed <- gData$map[rownames(gData$map) %in% colnames(markersRed), ]
      allFreq <- Matrix::colMeans(markersRed, na.rm = TRUE) / maxScore
      if (!useMAF) {
        MAF <- MAC / length(nonMiss) - 1e-5
      }
      ## Determine segregating markers. Exclude snps used as covariates.
      segMarkers <- which(allFreq >= MAF & allFreq <= (1 - MAF))
      ## Create data.frame for results.
      GWAResult <- data.frame(snp = rownames(mapRed), mapRed, pValue = NA,
                              LOD = NA, effect = NA, effectSe = NA, RLR2 = NA,
                              allFreq = allFreq, stringsAsFactors = FALSE)
      ## Define single column matrix with trait non missing values.
      y <- as(phEnvTr[which(phEnvTr$genotype %in% nonMiss), trait], "dgeMatrix")
      if (GLSMethod == 1) {
        ## Exclude snpCovariates from segregating markers.
        exclude <- computeExcludedMarkers(snpCov = snpCov,
                                          markersRed = markersRed,
                                          allFreq = allFreq)
        ## The following is based on the genotypes, not the replicates:
        X <- markersRed[nonMissRepId, setdiff(segMarkers, exclude)]
        if (length(covEnv) == 0) {
          Z <- NULL
        } else {
          ## Define covariate matrix Z.
          Z <- as(as.matrix(phEnvTr[which(phEnvTr$genotype %in% nonMiss),
                                    covEnv]), "dgeMatrix")
        }
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix, covs = Z)
        GWAResult[setdiff(segMarkers, exclude),
                  c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult
        ## Compute p-values and effects for snpCovariates using fastGLS.
        for (snpCovariate in snpCov) {
          GLSResultSnpCov <-
            fastGLS(y = y,
                    X = markersRed[nonMissRepId, snpCovariate, drop = FALSE],
                    Sigma = vcovMatrix,
                    covs = Z[, which(colnames(Z) != snpCovariate),
                             drop = FALSE],
                    nChunks = 1)
          GWAResult[snpCovariate, c("pValue", "effect", "effectSe", "RLR2")] <-
            GLSResultSnpCov
        }
      } else if (GLSMethod == 2) {
        ## Similar to GLSMethod 1 except using chromosome specific kinship
        ## matrices.
        for (chr in chrs) {
          mapRedChr <- mapRed[which(mapRed$chr == chr), ]
          markersRedChr <- markersRed[, which(colnames(markersRed) %in%
                                                rownames(mapRedChr)),
                                      drop = FALSE]
          allFreqChr <- Matrix::colMeans(markersRedChr, na.rm = TRUE) / maxScore
          ## Determine segregating markers. Exclude snps used as covariates.
          segMarkersChr <- which(allFreqChr >= MAF & allFreqChr <= (1 - MAF))
          ## Exclude snpCovariates from segregating markers.
          exclude <- computeExcludedMarkers(snpCov = snpCov,
                                            markersRed = markersRedChr,
                                            allFreq = allFreqChr)
          ## Remove excluded snps from segreg markers for current chromosome.
          segMarkersChr <- setdiff(intersect(segMarkersChr,
                                             which(mapRedChr$chr == chr)),
                                   exclude)
          X <- markersRedChr[nonMissRepId, segMarkersChr, drop = FALSE]
          if (length(covEnv) == 0) {
            Z <- NULL
          } else {
            ## Define covariate matrix Z.
            Z <- as(as.matrix(phEnvTr[which(phEnvTr$genotype %in% nonMiss),
                                      covEnv]),"dgeMatrix")
          }
          GLSResult <- fastGLS(y = y, X = X,
                               Sigma = vcovMatrix[[which(chrs == chr)]],
                               covs = Z)
          GWAResult[colnames(markersRedChr)[segMarkersChr],
                    c("pValue", "effect", "effectSe", "RLR2")] <-
            GLSResult
          ## Compute pvalues and effects for snpCovariates using fastGLS.
          for (snpCovariate in intersect(snpCov, colnames(markersRedChr))) {
            GLSResultSnpCov <-
              fastGLS(y = y, X = markersRed[nonMissRepId, snpCovariate,
                                            drop = FALSE],
                      Sigma = vcovMatrix[[which(chrs == chr)]],
                      covs = Z[, which(colnames(Z) != snpCovariate),
                               drop = FALSE], nChunks = 1)
            GWAResult[snpCovariate,
                      c("pValue", "effect", "effectSe", "RLR2")] <-
              GLSResultSnpCov
          }
        }
      }
      ## Effects should be for a single allele, not for 2
      if (maxScore == 1) {
        GWAResult$effect <- 0.5 * GWAResult$effect
      }
      ## Calculate the genomic inflation factor.
      GC <- genomicControlPValues(pVals = GWAResult$pValue,
                                  nObs = length(nonMiss), nCov = length(covEnv))
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
      ## When thrType is 1 or 3, determine the LOD threshold.
      if (thrType == 1) {
        ## Compute LOD threshold using Bonferroni correction.
        LODThr <- -log10(alpha / sum(!is.na(GWAResult$pValue)))
      } else if (thrType == 3) {
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
                     maxScore = maxScore, pheno = phEnvTr, trait = trait)
      GWATotEnv[[trait]] <- GWAResult
    } # end for (trait in traits)
    GWATot[[match(env, environments)]] <-
      as.data.frame(dplyr::bind_rows(GWATotEnv, .id = "trait"))
    signSnpTot[[match(env, environments)]] <-
      as.data.frame(dplyr::bind_rows(signSnpTotEnv, .id = "trait"))
    signSnpTot <- lapply(signSnpTot, FUN = function(x) {
      if (is.null(x) || nrow(x) == 0) NULL else x
    })
    varCompTot[[match(env, environments)]] <- varCompEnv
    LODThrTot[[match(env, environments)]] <- LODThrEnv
    inflationFactorTot[[match(env, environments)]] <- inflationFactorEnv
  } # end for (environment in environments)
  ## Collect info.
  GWASInfo <- list(call = match.call(),
                   remlAlgo = factor(remlAlgo, levels = c(1, 2),
                                     labels = c("EMMA", "Newton-Raphson")),
                   thrType = factor(thrType, levels = c(1, 2, 3),
                                    labels = c("bonferroni", "self-chosen",
                                               "smallest p-values")),
                   MAF = MAF,
                   GLSMethod =
                     factor(GLSMethod, levels = c(1, 2),
                            labels = c("single kinship matrix",
                                       "chromosome specific kinship matrices")),
                   varComp = varCompTot,
                   genomicControl = genomicControl,
                   inflationFactor = inflationFactorTot)
  return(createGWAS(GWAResult = GWATot,
                    signSnp = signSnpTot,
                    kin = if (GLSMethod == 1) {
                      K
                    } else {
                      KChr
                    },
                    thr = LODThrTot,
                    GWASInfo = GWASInfo))
}
