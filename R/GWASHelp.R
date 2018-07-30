#' Estimate variance components in single trait GWAS
#'
#' Helper function for estimating variance components in single trait GWAS.
#'
#' @keywords internal
estVarComp <- function(GLSMethod,
                       remlAlgo,
                       trait,
                       pheno,
                       covar,
                       K,
                       chrs,
                       KChr,
                       nonMiss,
                       nonMissRepId) {
  ## Estimate variance components.
  if (GLSMethod == 1) {
    if (isTRUE(all.equal(K, Matrix::Diagonal(nrow(K)), check.names = FALSE))) {
      ## Kinship matrix is computationally identical to identity matrix.
      vcovMatrix <- Matrix::Diagonal(nrow(pheno))
    }
  } else if (GLSMethod == 2) {
    varComp <- vcovMatrix <-
      setNames(vector(mode = "list", length = length(chrs)), paste("chr", chrs))
  }
  if (remlAlgo == 1) {
    ## emma algorithm takes covariates from gData.
    gDataEmma <-
      createGData(pheno = pheno[, c("genotype", trait)],
                  covar = if (is.null(covar)) {
                    NULL
                  } else {
                    as.data.frame(pheno[covar], row.names = pheno$genotype)
                  })
    if (GLSMethod == 1) {
      remlObj <- EMMA(gData = gDataEmma, trait = trait, environment = 1,
                      covar = covar, K = K)
      ## Extract varComp and vcovMatrix
      varComp <- remlObj$varComp
      vcovMatrix <- remlObj$vcovMatrix
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMiss, nonMiss]
        ## Compute variance components using chromosome specific kinship.
        remlObj <- EMMA(gData = gDataEmma, trait = trait,
                        environment = 1, covar = covar, K = K)
        ## Compute varcov matrix using var components.
        varComp[[which(chrs == chr)]] <- remlObj$varComp
        vcovMatrix[[which(chrs == chr)]] <- remlObj$vcovMatrix
      }
    }
  } else if (remlAlgo == 2) {
    if (!is.null(covar)) {
      ## Construct the formula for the fixed part of the model.
      ## Define formula for fixed part. ` needed to accommodate -
      ## in variable names.
      fixed <- as.formula(paste0(trait," ~ `",
                                 paste0(covar, collapse = "` + `"), "`"))
    } else {
      fixed <- as.formula(paste(trait, " ~ 1"))
    }
    if (GLSMethod == 1) {
      ## Fit model.
      modFit <- sommer::mmer2(fixed = fixed, data = pheno,
                              random = ~ g(genotype), G = list(genotype = K),
                              silent = TRUE, date.warning = FALSE)
      ## Compute varcov matrix using var components from model.
      vcMod <- modFit$var.comp
      modK <- K[nonMissRepId, nonMissRepId]
      varComp <- setNames(
        unlist(vcMod)[c(1, length(unlist(vcMod)))], c("Vg", "Ve"))
      vcovMatrix <- unlist(vcMod)[1] * modK +
        Matrix::Diagonal(n = nrow(modK),
                         x = unlist(vcMod)[length(unlist(vcMod))])
      if (any(eigen(vcovMatrix, symmetric = TRUE,
                    only.values = TRUE)$values <= 1e-8))
        vcovMatrix <- Matrix::nearPD(vcovMatrix)$mat
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMiss, nonMiss]
        ## Fit mmer2 model using chromosome specific kinship.
        modFit <- sommer::mmer2(fixed = fixed, data = pheno,
                                random = ~ g(genotype),
                                G = list(genotype = K), silent = TRUE,
                                date.warning = FALSE)
        ## Compute varcov matrix using var components from model.
        vcMod <- modFit$var.comp
        modK <- K[nonMissRepId, nonMissRepId]
        varComp[[which(chrs == chr)]] <- setNames(
          unlist(vcMod)[c(1, length(unlist(vcMod)))], c("Vg", "Ve"))
        vcovMatrix[[which(chrs == chr)]] <- unlist(vcMod)[1] * modK +
          unlist(vcMod)[length(unlist(vcMod))] *
          Matrix::Diagonal(n = nrow(modK))
      }
      vcovMatrix <- lapply(vcovMatrix, FUN = function(vc) {
        if (any(eigen(vc, symmetric = TRUE,
                      only.values = TRUE)$values <= 1e-8)) {
          Matrix::nearPD(vc)$mat
        } else {
          vc
        }
      })
    }
  }
  return(list(varComp = varComp, vcovMatrix = vcovMatrix))
}

#' Select markers to be excluded from GWAS scan.
#'
#' Helper function for selecting markers to be excluded from GWAS scan.
#' Markers are excluded if they are identical to any of the snpCovariates
#' (including the snpCovariates themselves).
#'
#' @param snpCov A character vector of snpCovariates.
#' @param markers A matrix with marker information.
#' @param allFreq A numerical vector of allele frequencies of the markers in
#' \code{markers}. This could be computed from markers as well but it is
#' needed in the general algorithm so to not redo things unnecessarily it is
#' not redone here.
#'
#' @return A numerical vector of markers to be exluded from the GWAS scan.
#'
#' @keywords internal
exclMarkers <- function(snpCov,
                        markers,
                        allFreq) {
  exclude <- integer()
  if (any(snpCov %in% colnames(markers))) {
    snpCovNumbers <- which(colnames(markers) %in% snpCov)
    for (snp in snpCovNumbers) {
      ## Rough selection based on allele frequency. Done for speed.
      candidates <- which(allFreq == allFreq[snp])
      ## Exclude all snps that are identical to snps in snpCovariates.
      snpInfo <- as.numeric(markers[, snp])
      exclude <- union(exclude,
                       candidates[apply(X = markers[, candidates, drop = FALSE],
                                        MARGIN = 2, FUN = function(x) {
                                          identical(as.numeric(x), snpInfo)
                                        })])
    }
  }
  return(exclude)
}

#' Correction of p-values based on genomic inflation
#'
#' Correction of p-values based on the genomic inflation factor, as in Devlin
#' and Roeder (1999). It is assumed that the p-values come from an F-test with
#' df1 = 1 and df2 = nObs - nCov - 2.
#'
#' @param pVals A numeric vector of p-values between 0 and 1; may contain NA's.
#' @param nObs An integer > 0 indicating the number of individuals.
#' @param nCov An integer > 0 indicating the number of covariables.
#'
#' @return A list with two components:
#' \itemize{
#' \item{\code{pValues} a vector of p-values corrected by the genomic inflation
#' factor, with the same NA's as the input}.
#' \item{\code{inflation} the inflation factor}.
#' }
#'
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association
#' studies. Biometrics, December 1999, Vol. 55(4), p. 997-1004.
#'
#' @keywords internal
genCtrlPVals <- function(pVals,
                         nObs,
                         nCov = 0) {
  ## Check input.
  if (missing(pVals) || !is.numeric(pVals) || any(pVals < 0, na.rm = TRUE) ||
      any(pVals > 1, na.rm = TRUE)) {
    stop("pVals should be a numeric vector with values between 0 and 1.\n")
  }
  if (missing(nObs) || length(nObs) > 1 || !is.numeric(nObs) ||
      nObs != round(nObs) || nObs < 1) {
    stop("nObs should be a single positive integer.\n")
  }
  if (length(nCov) > 1 || !is.numeric(nCov) || nCov != round(nCov) ||
      nCov < 0) {
    stop("nCov should be a single non negative integer.\n")
  }
  ## Compute degree of freedom.
  df2 <- nObs - nCov - 2
  pValsNew <- pVals
  ## Compute F-values from input p-values.
  fVals <- qf(p = na.omit(pVals), df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute inflation factor as in Devlin and Roeder.
  inflation <- median(fVals, na.rm = TRUE) / qf(p = 0.5, df1 = 1, df2 = df2,
                                                lower.tail = FALSE)
  ## Compute new F-values and p-values.
  fValsNew <- fVals / inflation
  pValsNew[!is.na(pVals)] <- pf(q = fValsNew, df1 = 1, df2 = df2,
                                lower.tail = FALSE)
  return(list(pValues = pValsNew, inflation = inflation))
}

#' @keywords internal
extrSignSnps <- function(GWAResult,
                         LODThr,
                         sizeInclRegion,
                         minR2,
                         map,
                         markers,
                         maxScore,
                         pheno,
                         trait) {
  signSnpNr <- which(!is.na(GWAResult$LOD) & GWAResult$LOD >= LODThr)
  if (length(signSnpNr) > 0) {
    if (sizeInclRegion > 0) {
      snpSelection <-
        unlist(sapply(X = signSnpNr, FUN = getSNPsInRegionSufLD,
                      ## Create new minimal gData object to match map and
                      ## markers used for SNP selection.
                      gData = createGData(map = map, geno = markers),
                      regionSize = sizeInclRegion, minR2 = minR2))
      snpSelection <- sort(union(snpSelection, signSnpNr))
      snpStatus <- rep(paste("within", sizeInclRegion / 1000,
                             "kb of a significant snp"),
                       length(snpSelection))
      snpStatus[snpSelection %in% signSnpNr] <- "significant snp"
    } else {
      snpSelection <- signSnpNr
      snpStatus <- rep("significant snp", length(signSnpNr))
    }
    effect <- GWAResult$effect[snpSelection]
    ## Compute variance of marker scores, based on genotypes for which
    ## phenotypic data is available. For inbreeders, this depends on
    ## maxScore. It is therefore scaled to marker scores 0, 1 (or 0, 0.5,
    ## 1 if there are heterozygotes)
    snpVar <- 4 * effect ^ 2 / maxScore ^ 2 *
      apply(X = markers[, snpSelection, drop = FALSE], MARGIN = 2, FUN = var)
    propSnpVar <- snpVar / as.numeric(var(pheno[trait]))
    ## Create data.frame with significant snps.
    GWAResultSel <- GWAResult[snpSelection, ]
    signSnp <- data.frame(GWAResult[snpSelection, ],
                          snpStatus = as.factor(snpStatus),
                          propSnpVar = propSnpVar, stringsAsFactors = FALSE)
  } else {
    signSnp <- data.frame()
  }
  return(signSnp)
}

#' get the SNPs close to a given SNP with sufficient LD
#'
#' \code{getSNPsInRegionSufLD} extracts the SNPs from a map file that are
#' within a given distance of a reference SNP (on either side). Only those SNPs
#' that are in sufficient linkage disequilibrium (LD) with the reference SNP
#' are returned.
#'
#' @param gData An object of class gData with at least the map and markers
#' included.
#' @param snp An integer indicating the index of the reference SNP within
#' the map.
#' @param regionSize A numerical value indicating the size of the region on
#' the chromosome in which to look for SNPs.
#' @param minR2 A numerical value between 0 and 1 indicating the minimum
#' LD (in terms of r^2) that the SNPs should have with the reference SNP.
#'
#' @return An integer vector with indices of the SNPs that are within the
#' given \code{regionSize} and have a minimum LD with the reference SNP.
#'
#' @keywords internal
getSNPsInRegionSufLD <- function(gData,
                                 snp,
                                 regionSize = 5000,
                                 minR2 = 0.5) {
  ## Check input.
  if (missing(gData) || !is.gData(gData) || is.null(gData$map) ||
      is.null(gData$markers)) {
    stop(paste("gData should be a valid gData object containing at least",
               "map and markers.\n"))
  }
  if (missing(snp) || length(snp) > 1 || !is.numeric(snp) ||
      snp != round(snp) || !snp %in% 1:nrow(gData$map)) {
    stop(paste("snp should be a single integer indicating a row in",
               "the map in gData.\n"))
  }
  if (length(regionSize) > 1 || !is.numeric(regionSize) || regionSize < 0) {
    stop("regionSize should be a single positive numerical value.\n")
  }
  if (length(minR2) > 1 || !is.numeric(minR2) || minR2 < 0 || minR2 > 1) {
    stop("minR2 should be a single numerical value between 0 and 1.")
  }
  ## Get candidate SNPs based on position.
  crit1 <- abs(gData$map$pos[snp] - gData$map$pos) <= regionSize
  crit2 <- gData$map$chr == gData$map$chr[snp]
  candidateSnps <- setdiff(which(crit1 & crit2), snp)
  ## Compute R2 for candidate SNPs.
  R2 <- suppressWarnings(cor(as.matrix(gData$markers[, candidateSnps]),
                             gData$markers[, snp]) ^ 2)
  ## Select SNPs based on R2.
  candidateSnpsNames <- names(which(R2[, 1] > minR2))
  return(which(rownames(gData$map) %in% candidateSnpsNames))
}
