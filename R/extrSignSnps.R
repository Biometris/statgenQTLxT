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

