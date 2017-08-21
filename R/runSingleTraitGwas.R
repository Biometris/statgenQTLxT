#' perform single-trait GWAS
#'
#' \code{runSingleTraitGwas} performs a single-trait Genome Wide Association Study (GWAS) on phenotypic and
#' genotypic data contained in a \code{gData} object. A covariance matrix is computed using the EMMA algorithm
#' (Kang et al., 2008) or the Newton-Raphson algorithm (Tunnicliffe, 1989) of the \code{sommer}
#' package. Then a Generalized Least Squares (GLS) method is used for estimating the marker effects and
#' corresponding p-values. This is done using either one kinship matrix for all chromosomes or different
#' chromosome specific kinship matrices for each chromosome. Significant SNPs are selected based on a user
#' defined threshold.
#'
#' @param gData an object of class \code{gData} containing at least \code{map}, \code{markers} and
#' \code{pheno}.
#' @param traits a vector of traits on which to run GWAS. These can be either numeric indices or
#' character names of columns in \code{pheno}. If \code{NULL} GWAS is run on all traits.
#' @param environments a vector of environments on which to run GWAS. These can be either numeric indices or
#' character names of list items in \code{pheno}. If \code{NULL} GWAS is run for all environments.
#' @param covar an optional vector of covariates taken into account when running GWAS. These can be either
#' numeric indices or character names of columns in \code{covar} in \code{gData}. If \code{NULL} no
#' covariates are used.
#' @param snpCovariates an optional character vector of snps to be included as covariates.
#' @param K an optional kinship matrix or list of kinship matrices. If \code{GLSMethod} = 1 then one
#' matrix should be provided, if \code{GLSMethod} = 2 a list of chromosome specific matrices of lenght
#' equal to the number of chromosomes in \code{map} in \code{gData}. If \code{NULL} then matrix
#' \code{kinship} in \code{gData} is used. If both \code{K} is provided and \code{gData} contains a
#' matrix \code{kinship} then \code{K} is used.
#' @param kinshipMethod an optional character indicating the method used for calculating the kinship
#' matrix(ces). Currently "astle" (Astle and Balding, 2009), "GRM", "IBS" and "vanRaden" (VanRaden, 2008)
#' are supported. If a kinship matrix is supplied either in \code{gData} or in parameter \code{K}
#' \code{kinshipMatrix} is ignored.
#' @param remlAlgo an integer indicating the algorithm used to estimate the variance components
#' \enumerate{
#' \item{EMMA}
#' \item{Newton-Raphson using sommer package}
#' }
#' @param GLSMethod an integer indicating the method used to estimate the marker effects.
#' \enumerate{
#' \item{using one kinship matrix.}
#' \item{using chromosome specific kinship matrices; similar to the approach of Rincent et al.}
#' }
#' @param sizeInclRegion an integer. Should the results for SNPs close to significant SNPs
#' be included? If so, the size of the region in centimorgan or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included in the region close
#' to significant SNPs to only those SNPs that are in sufficient Linkage Disequilibruim (LD) with
#' the significant snp, where LD is measured in terms of r^2. If for example
#' \code{sizeInclRegion} = 200000 and \code{minR2} = 0.5, then for every significant SNP
#' also those SNPs whose LD (r^2) with the significant SNP is at least 0.5 AND which are
#' at most 200kb away from this significant snp are included. Ignored if \code{sizeInclRegion} = 0.
#' @param useMAF should the minor allele frequency be used for selecting SNPs for the analysis.
#' If \code{FALSE} the minor allele count is used instead.
#' @param MAF A numerical value between 0 and 1. SNPs with minor allele frequency below
#' this value are not taken into account for the analysis, i.e. p-values and effect sizes are
#' put to missing (NA). Ignored if \code{useMAF} is \code{FALSE}.
#' @param MAC A numerical value. SNPs with minor allele count below this value are not
#' taken into account for the analysis, i.e. p-values and effect sizes are put to missing (NA).
#' Ignored if \code{useMAF} is \code{TRUE}.
#' @param genomicControl should genomic control correction as in Devlin and Roeder (1999) be applied?
#' @param thrType an integer indicating the type of threshold used for the selection of candidate
#' loci.
#' \enumerate{
#' \item{Bonferroni; a LOD-threshold of \eqn{-log10(alpha/p)}, where p is the number of markers and alpha
#' can be specified in \code{alpha}.}
#' \item{a self-chosen LOD-threshold, specied in \code{LODThr}.}
#' \item{the LOD-threshold is chosen such the the SNPs with the \code{nSnpLOD} smallest p-values are selected.
#' \code{nSnpLOD} can be specified.}
#' }
#' @param alpha a numerical value used for calculating the LOD-threshold for \code{thrType} = 1.
#' @param LODThr a numerical value used as a LOD-threshold when \code{thrType} = 2.
#' @param nSnpLOD a numerical value indicating the number of SNPs with the smallest p-values that
#' are selected when \code{thrType} = 3.
#'
#' @return an object of \code{\func{class GWAS}}
#'
#' @references Kang et al. (2008) Efficient Control of Population Structure in Model Organism
#' Association Mapping. Genetics, March 2008, Vol. 178, no. 3, p. 1709-1723
#' @references Tunnicliffe W. (1989) On the use of marginal likelihood in time series model estimation.
#' JRSS, Vol.51(1), p.15-27.
#' @references Rincent et al. (2014) Recovering power in association mapping panels with variable
#' levels of linkage disequilibrium. Genetics, May 2014. Vol. 197. p. 375–387.
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association studies. Biometrics,
#' December 1999, Vol. 55(4), p. 997-1004.
#' @references Segura et al. (2012) An efficient multi-locus mixed-model approach for genome-wide
#' association studies in structured populations. Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association mapping.
#' Heredity, February 2010, Vol. 105, p. 333–340.
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic relatedness in genetic
#' association studies, Stat. Sci., November 2009, Vol. 24, no. 4, p. 451–471.L
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic predictions. J Dairy Sci,
#' November 2008, Vol. 91 p. 4414–4423.
#'
#' @import stats

runSingleTraitGwas <- function (gData,
  traits = NULL,
  environments = NULL,
  covar = NULL,
  snpCovariates = NULL,
  K = NULL,
  kinshipMethod = "astle",
  remlAlgo = 1,
  GLSMethod = 1,
  sizeInclRegion = 0,
  minR2,
  useMAF = TRUE,
  MAF = 0.05,
  MAC = 10,
  genomicControl = FALSE,
  thrType = 1,
  alpha = 0.05 ,
  LODThr = 4,
  nSnpLOD = 10) {
  ## Checks.
  if(missing(gData) || !is.gData(gData) || is.null(gData$map) || is.null(gData$markers) ||
      is.null(gData$pheno))
    stop("gData should be a valid gData object with at least map, markers and pheno included.\n")
  if (!is.null(environments) && !is.numeric(environments) && !is.character(environments))
    stop("environments should be a numeric or character vector.\n")
  if ((is.character(environments) && !all(environments %in% names(gData$pheno))) ||
      (is.numeric(environments) && any(environments > length(gData$pheno))))
    stop("environments should be list items in pheno.\n")
  if (!is.null(traits) && !is.numeric(traits) && !is.character(traits))
    stop("traits should be a numeric or character vector.\n")
  for (environment in environments) {
    if ((is.character(traits) && !all(traits %in% colnames(gData$pheno[[environment]]))) ||
        (is.numeric(traits) && (any(traits == 1) || any(traits > ncol(gData$pheno[[environment]])))))
      stop("traits should be columns in pheno.\n")
  }
  ## If environments is null set environments to all environments in pheno.
  if (is.null(environments)) environments <- 1:length(gData$pheno)
  if (!is.null(covar) && !is.numeric(covar) && !is.character(covar))
    stop("covar should be a numeric or character vector.\n")
  if ((is.character(covar) && !all(covar %in% colnames(gData$covar))) ||
      (is.numeric(covar) && any(covar > ncol(gData$covar))))
    stop("covar should be columns in covar.\n")
  ## If covar is given as numeric convert to character.
  if (is.numeric(covar)) covar <- colnames(gData$covar)[covar]
  if (!is.null(snpCovariates) && !all(snpCovariates %in% colnames(gData$markers)))
    stop("All snpCovariates should be in markers.\n")
  if (is.null(remlAlgo) || length(remlAlgo) > 1 || !is.numeric(remlAlgo))
    stop("remlAlgo should be a single numeric value.\n")
  if (is.null(GLSMethod) || length(GLSMethod) > 1 || !is.numeric(GLSMethod))
    stop("GLSMethod should be a single numeric value.\n")
  if (GLSMethod == 1 && !is.null(K) && !is.matrix(K))
    stop("K should be a matrix.\n")
  if (GLSMethod == 2 && !is.null(K) && (!is.list(K) ||
      !all(sapply(K, FUN = is.matrix)) || length(K) != length(unique(gData$map$chr))))
    stop("K should be a list of matrices of length equal to the number of chromosomes
      in the map.\n")
  if ((GLSMethod == 1 && is.null(gData$kinship) && is.null(K)) ||
      GLSMethod == 2 && is.null(K))
    kinshipMethod <- match.arg(kinshipMethod, choices = c("astle", "GRM", "IBS", "vanRaden"))
  if (is.null(sizeInclRegion) || length(sizeInclRegion) > 1 || !is.numeric(sizeInclRegion) ||
      round(sizeInclRegion) != sizeInclRegion)
    stop("sizeInclRegion should be a single integer\n")
  if (sizeInclRegion > 0) {
    if (missing(minR2) || length(minR2) > 1 || !is.numeric(minR2) || minR2 < 0 || minR2 > 1)
      stop("minR2 should be a single numerical value between 0 and 1.\n")
  }
  if (useMAF) {
    if (missing(MAF) || length(MAF) > 1 || !is.numeric(MAF) || MAF < 0 || MAF > 1)
      stop("MAF should be a single numerical value between 0 and 1.\n")
    if (MAF == 0) {MAF <- 1e-6}
  } else {
    if (missing(MAC) || length(MAC) > 1 || !is.numeric(MAC))
      stop("MAF should be a single numerical value.\n")
    if (MAC == 0) {MAC <- 1}
  }
  if (is.null(thrType) || length(thrType) > 1 || !is.numeric(thrType))
    stop("thrType should be a single numeric value.\n")
  if (thrType == 1) {
    if (is.null(alpha) || length(alpha) > 1 || !is.numeric(alpha))
      stop("alpha should be a single numerical value.\n")
  } else if (thrType == 2) {
    if (is.null(LODThr) || length(LODThr) > 1 || !is.numeric(LODThr))
      stop("LODThr should be a single numerical value.\n")
  } else if (thrType == 3) {
    if (is.null(nSnpLOD) || length(nSnpLOD) > 1 || !is.numeric(nSnpLOD))
      stop("nSnpLOD should be a single numerical value.\n")
  }
  if (!(remlAlgo %in% 1:2)) {
    remlAlgo <- 2
    warning("Invalid value for remlAlgo. remlAlgo set to 2.\n")
  }
  if (!(GLSMethod %in% 1:2)) {
    GLSMethod <- 1
    warning("Invalid value for GLSMethod. GLSMethod set to 1.\n")
  }
  if (!(thrType %in% 1:3)) {
    thrType <- 1
    warning("Invalid value for thrType. thrType set to 1.\n")
  }
  ## Compute kinship matrix.
  if (GLSMethod == 1 && is.null(gData$kinship) && is.null(K))
      K <- do.call(kinshipMethod, list(X = gData$markers))
  ## Compute kinship matrices per chromosome. Only needs to be done once.
  if (GLSMethod == 2) {
    chrs <- unique(gData$map$chr[rownames(gData$map) %in% colnames(gData$markers)])
    if (!is.null(K)) {
      ## K is supplied. Set KChr to K.
      KChr <- K
    } else {
      ## Create list of zero matrices.
      KChr <- replicate(length(chrs),
        matrix(data = 0, nrow = nrow(gData$markers), ncol = nrow(gData$markers)), simplify = FALSE)
      for (chr in chrs) {
        ## Extract markers for current chromosome.
        chrMrk <- which(colnames(gData$markers) %in% rownames(gData$map[gData$map$chr == chr, ]))
        ## Compute kinship for current chromosome only.
        K <- do.call(kinshipMethod, list(X = gData$markers[, chrMrk], denominator = ncol(gData$markers)))
        ## Add computed kinship to all other matrices in KChr.
        for (i in setdiff(1:length(chrs), which(chr == chrs))) KChr[[i]] <- KChr[[i]] + K
      }
    }
  }
  ## Compute max value in markers
  maxScore <- max(gData$markers, na.rm = TRUE)
  # allele frequencies based on all genotypes (trait-independent)
  allFreqTot <- colMeans(gData$markers, na.rm = TRUE)
  if (maxScore == 2) {
    allFreqTot<- allFreqTot / 2
  }
  ## Define data.frames for total results.
  GWATot <- signSnpTot <- vector(mode = "list", length = length(environments))
  for (environment in environments) {
    ## If traits is given as numeric convert to character.
    if (is.numeric(traits)) traits <- colnames(gData$pheno[[environment]])[traits]
    ## If no traits supplied extract them from pheno data.
    if (is.null(traits)) traits <- colnames(gData$pheno[[environment]])[-1]
    ## Add covariates to pheno data.
    if (is.null(covar)) {
      phenoEnvir <- gData$pheno[[environment]]
      covarEnvir <- NULL
    } else {
      ## Append covariates to pheno data. Merge to remove values from pheno that are missing in covar.
      phenoEnvir <- merge(gData$pheno[[environment]], gData$covar[covar], by.x = "genotype", by.y = "row.names")
      ## Remove rows from phenoEnvir with missing covar check if there are missing values.
      phenoEnvir <- phenoEnvir[!is.na(phenoEnvir[covar]), ]
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
    GWATotEnvir <- signSnpTotEnvir <- NULL
    ## Perform GWAS for all traits.
    for (trait in traits) {
      ## Select relevant columns only.
      phenoEnvirTrait <- phenoEnvir[!is.na(phenoEnvir[trait]), c("genotype", trait, covarEnvir)]
      ## Select genotypes where trait is not missing.
      nonMissing <- unique(phenoEnvirTrait$genotype)
      if (GLSMethod == 1) {
        if (is.null(K)) {
          kinshipRed <- gData$kinship[nonMissing, nonMissing]
        } else {
          kinshipRed <- K[nonMissing, nonMissing]
        }
      }
      nonMissingRepId <- phenoEnvirTrait$genotype
      ## Estimate variance components.
      if (GLSMethod == 1) {
        if (!isTRUE(all.equal(kinshipRed, diag(nrow(kinshipRed)), check.names = FALSE))) {
          if (remlAlgo == 1) {
            ## emma algorithm takes covariates from gData.
            gDataEmma <- createGData(pheno = gData$pheno,
              covar = as.data.frame(phenoEnvir[covarEnvir], row.names = phenoEnvir$genotype))
            remlObj <- runEmma(gData = gDataEmma, trait = trait, environment = environment,
              covar = covarEnvir, K = kinshipRed)
            ## Compute varcov matrix using var components.
            varComp <- remlObj[[1]]
            vcovMatrix <- remlObj[[1]][1] * remlObj[[2]] +
              remlObj[[1]][2] * diag(nrow(remlObj[[2]]))
          } else if (remlAlgo == 2) {
            ## Construct the formula for the fixed part of the model.
            if (!is.null(covarEnvir)) {
              ## Define formula for fixed part. ` needed to accommodate - in variable names.
              fixed <- as.formula(paste0(trait," ~ `", paste0(covarEnvir, collapse = "` + `"), "`"))
            } else {
              fixed <- as.formula(paste(trait, " ~ 1"))
            }
            ## Fit mmer2 model.
            sommerFit <- sommer::mmer2(fixed = fixed, data = phenoEnvirTrait,
              random = ~ sommer::g(genotype), G = list(genotype = kinshipRed), silent = TRUE)
            ## Compute varcov matrix using var components from model.
            varComp <- sommerFit$var.comp[c(1, nrow(sommerFit$var.comp)), 1]
            vcovMatrix <- sommerFit$var.comp[1, 1] * kinshipRed +
              diag(sommerFit$var.comp[nrow(sommerFit$var.comp), 1], nrow = nrow(kinshipRed))
          }
        } else {
          ## Kinship matrix is computationally identical to identity matrix.
          vcovMatrix <- diag(nrow(phenoEnvirTrait))
        }
      } else if (GLSMethod == 2) {
        varComp <- vcovMatrix <- vector(mode = "list", length = length(chrs))
        names(varComp) <- paste("chr", chrs)
        ## emma algorithm takes covariates from gData.
        if (remlAlgo == 1) {
          gDataEmma <- createGData(pheno = gData$pheno,
            covar = as.data.frame(phenoEnvir[covarEnvir], row.names = phenoEnvir$genotype))
          for (chr in chrs) {
            ## Get chromosome specific kinship.
            KinshipRedChr <- KChr[[which(chrs == chr)]][nonMissing, nonMissing]
            ## Compute variance components using chromosome specific kinship.
            remlObj <- runEmma(gData = gDataEmma, trait = trait, environment = environment,
              covar = covarEnvir, K = KinshipRedChr)
            ## Compute varcov matrix using var components.
            varComp[[which(chrs == chr)]] <- remlObj[[1]]
            vcovMatrix[[which(chrs == chr)]] <- remlObj[[1]][1] * remlObj[[2]] +
              remlObj[[1]][2] * diag(nrow(remlObj[[2]]))
          }
        } else if (remlAlgo == 2) {
          if (!is.null(covarEnvir)) {
            ## Define formula for fixed part. ` needed to accommodate - in variable names.
            fixed <- as.formula(paste0(trait," ~ `", paste0(covarEnvir, collapse = "` + `"), "`"))
          } else {
            fixed <- as.formula(paste(trait, " ~ 1"))
          }
          for (chr in chrs) {
            ## Get chromosome specific kinship.
            KinshipRedChr <- KChr[[which(chrs == chr)]][nonMissing, nonMissing]
            ## Fit mmer2 model using chromosome specific kinship.
            sommerFit <- sommer::mmer2(fixed = fixed, data = phenoEnvirTrait,
              random = ~ sommer::g(genotype), G = list(genotype = KinshipRedChr), silent = TRUE)
            ## Compute varcov matrix using var components from model.
            varComp[[which(chrs == chr)]] <- sommerFit$var.comp[c(1, nrow(sommerFit$var.comp)), 1]
            vcovMatrix[[which(chrs == chr)]] <- sommerFit$var.comp[1, 1] * KinshipRedChr +
              diag(sommerFit$var.comp[nrow(sommerFit$var.comp), 1], nrow = nrow(KinshipRedChr))
          }
        }
      }
      ## Compute allele frequencies based on genotypes for which phenotypic data is available.
      markersRed <- gData$markers[nonMissing, ]
      mapRed <- gData$map[rownames(gData$map) %in% colnames(markersRed), ]
      allFreq <- colMeans(markersRed)
      if (maxScore == 2) {
        allFreq <- allFreq / 2
      }
      if (!useMAF) {MAF <- MAC / length(nonMissing) - 1e-5}
      ## Determine segregating markers. Exclude snps used as covariates.
      segMarkers <- which(allFreq >= maxScore * MAF & allFreq <= maxScore * (1 - MAF))
      ## Exclude snpCovariates from segragating markers.
      if (!is.null(snpCovariates)) {
        snpCovariateNumbers <- which(colnames(markersRed) %in% snpCovariates)
        for (snp in snpCovariateNumbers) {
          ## Rough selection based on allele frequency. Done for speed.
          candidates <- which(allFreq == allFreq[snp])
          ## Exclude all snps that are identical to snps in snpCovariates.
          snpInfo <- markersRed[, snp]
          exclude <- candidates[apply(markersRed[, candidates], MARGIN = 2,
            FUN = function(x) {identical(as.numeric(x), as.numeric(snpInfo))})]
        }
        ## Remove excluded snps from segregating markers.
        segMarkers <- setdiff(segMarkers, exclude)
      }
      ## Define data.frame for results.
      GWAResult <- data.frame(snp = rownames(mapRed),
        mapRed,
        pValue = NA,
        LOD = NA,
        effect = NA,
        effectSe = NA,
        RLR2 = NA,
        allFreq = allFreq,
        stringsAsFactors = FALSE)
      y <- phenoEnvirTrait[which(phenoEnvirTrait$genotype %in% nonMissing), trait]
      if (GLSMethod == 1) {
        ## The following is based on the genotypes, not the replicates:
        X <- markersRed[nonMissingRepId, segMarkers]
        if (length(covarEnvir) == 0) {
          ## Compute pvalues and effects using fastGLS.
          GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix)
        } else {
          ## Define covariate matrix Z.
          Z <- as.matrix(phenoEnvirTrait[which(phenoEnvirTrait$genotype %in% nonMissing), covarEnvir])
          ## Compute pvalues and effects using fastGLS.
          GLSResult <- fastGLSCov(y = y, X = X, Sigma = vcovMatrix, covs = Z)
        }
        GWAResult[segMarkers, c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult
      } else if (GLSMethod == 2) {
        ## Same as for GLSMethod 1 but then per chromosome using the chromosome specific
        ## kinship calculated before.
        for (chr in chrs) {
          segMarkersChr <- intersect(segMarkers, which(mapRed$chr == chr))
          X <- markersRed[nonMissingRepId, segMarkersChr]
          if (length(covarEnvir) == 0) {
            GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix[[which(chrs == chr)]])
          } else {
            Z <- as.matrix(phenoEnvirTrait[which(phenoEnvirTrait$genotype %in% nonMissing), covarEnvir])
            GLSResult <- fastGLSCov(y = y, X = X, Sigma = vcovMatrix[[which(chrs == chr)]], covs = Z)
          }
          GWAResult[segMarkersChr, c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult
        }
      }
      ## Effects should be for a single allele, not for 2
      if (maxScore == 1) {
        GWAResult$effect <- 0.5 * GWAResult$effect
      }
      ## Calculate the genomic inflation factor and rescale p-values.
      if (genomicControl) {
        GC <- genomicControlPValues(pVals = GWAResult$pValue,
          nObs = length(nonMissing),
          nCov = length(covarEnvir))
        GWAResult$pValue <- GC[[1]]
      }
      ## Compute LOD score.
      GWAResult$LOD <- -log10(GWAResult$pValue)
      ## Add gene information if available.
      if ("genes" %in% names(gData)) {
        GWAResult <- cbind(GWAResult, gene1 = gData$genes$gene1, gene2 = gData$genes$gene2)
      }
      ## When thrType is 1, 3 or 4, determine the LOD threshold.
      if (thrType == 1) {
        ## Compute LOD threshold using Bonferroni correction.
        nEff <- sum(!is.na(GWAResult$pValue))
        LODThr <- -log10(alpha / nEff)
      } else if (thrType == 3) {
        ## Compute LOD threshold by computing the 10log of the nSnpLOD item of ordered p values.
        LODThr <- -log10(sort(na.omit(GWAResult$pValue))[nSnpLOD])
      }
      ## Select the SNPs whose LOD-scores is above the threshold
      signSnp <- which(!is.na(GWAResult$pValue) & -log10(GWAResult$pValue) >= LODThr)
      if (length(signSnp) > 0) {
        if (sizeInclRegion > 0) {
          snpSelection <- unlist(sapply(signSnp,
            FUN = getSNPsInRegionSufLD,
            ## Create new minimal gData object to match map and markers used for SNP selection.
            gData = createGData(map = mapRed, geno = markersRed),
            regionSize = sizeInclRegion,
            minR2 = minR2))
          snpSelection <- sort(union(snpSelection, signSnp))
          snpStatus <- rep(paste("within", sizeInclRegion/1000 , "kb of a significant snp"),
            length(snpSelection))
          snpStatus[snpSelection %in% signSnp] <- "significant snp"
        } else {
          snpSelection <- signSnp
          snpStatus <- rep("significant snp", length(signSnp))
        }
        if (GLSMethod %in% 1:2) {
          effects <- GWAResult$effect[snpSelection]
          ## Compute variance of marker scores, based on genotypes for which phenotypic data is available.
          ## for inbreeders, this depends on maxScore. It is therefore scaled to marker scores 0, 1
          ## (or 0, 0.5, 1 if there are heterozygotes)
          snpVar <- 4 * effects ^ 2 *
            apply(as.matrix(markersRed[, snpSelection]), MARGIN = 2, FUN = var) / maxScore ^ 2
          propSnpVar <- snpVar / as.numeric(var(phenoEnvirTrait[trait]))
        }
        ## Create data.frame with significant snps.
        signSnp <- data.frame(snp = colnames(gData$markers)[snpSelection],
          gData$map[snpSelection, ],
          pValue = GWAResult$pValue[snpSelection],
          LOD = GWAResult$LOD[snpSelection],
          snpStatus = as.factor(snpStatus),
          allFreq = allFreq[snpSelection],
          effects,
          RLR2 = GWAResult$RLR2[snpSelection],
          propSnpVar = propSnpVar,
          stringsAsFactors = FALSE)
        signSnpTotEnvir <- rbind(signSnpTotEnvir, data.frame(trait = trait, signSnp, stringsAsFactors = FALSE))
      }
      GWATotEnvir <- rbind(GWATotEnvir, data.frame(trait = trait, GWAResult, stringsAsFactors = FALSE))
    } # end for (trait in traits)
    GWATot[[match(environment, environments)]] <- GWATotEnvir
    signSnpTot[[match(environment, environments)]] <- signSnpTotEnvir
  } # end for (environment in environments)
  names(GWATot) <- names(signSnpTot) <- names(gData$pheno[environments])
  ## Collect info
  GWASInfo <- list(call = match.call(),
    GLSMethod = factor(GLSMethod, levels = c(1, 2), labels = c("EMMA", "Newton-Raphson")),
    thrType = factor(thrType, levels = c(1, 2, 3),
      labels = c("bonferroni", "self-chosen", "smallest p-values")),
    MAF = MAF,
    varComp = varComp,
    genomicControl = genomicControl)
  if (genomicControl) GWASInfo$inflationFactor <- GC[[2]]
  return(createGWAS(GWAResult = GWATot,
    signSnp = signSnpTot,
    kin = if (GLSMethod == 1) K else KChr,
    thr = LODThr,
    GWASInfo = GWASInfo))
}
