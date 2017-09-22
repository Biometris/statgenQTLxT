#' get the SNPs close to a given SNP with sufficient LD
#'
#' \code{getSNPsInRegionSufLD} extracts the SNPs from a map file that are within a given distance
#' of a reference SNP (on either side). Only those SNPs that are in sufficient linkage disequilibrium
#' (LD) with the reference SNP are returned.
#'
#' @param gData an object of class gData with at least the map and markers included.
#' @param snp an integer indicating the index of the reference SNP within the map.
#' @param regionSize a numeric value indicating the size of the region on the chromosome in which to look
#' for SNPs.
#' @param minR2 a numeric value between 0 and 1 indicating the minimum LD (in terms of r^2) that the SNPs
#' should have with the reference SNP.
#'
#' @return an integer vector with indices of the SNPs that are within the given \code{regionSize} and
#' have a minimum LD with the reference SNP.
#'
#' @import stats
#'
#' @keywords internal

getSNPsInRegionSufLD <- function(gData,
  snp,
  regionSize = 5000,
  minR2 = 0.5) {
  ## Check input.
  if (missing(gData) || !is.gData(gData) || is.null(gData$map) || is.null(gData$markers))
    stop("gData should be a valid gData object containing at least map and markers.\n")
  if (missing(snp) || length(snp) > 1 || !is.numeric(snp) || snp != round(snp) ||
      !snp %in% 1:nrow(gData$map))
    stop("snp should be a single integer indicating a row in the map in gData.\n")
  if (length(regionSize) > 1 || !is.numeric(regionSize) || regionSize < 0)
    stop("regionSize should be a single positive numeric value.\n")
  if (length(minR2) > 1 || !is.numeric(minR2) || minR2 < 0 || minR2 > 1)
    stop("minR2 should be a single numeric value between 0 and 1.")
  ## Get candidate SNPs based on position.
  crit1 <- abs(gData$map$pos[snp] - gData$map$pos) <= regionSize
  crit2 <- gData$map$chr == gData$map$chr[snp]
  candidateSnps <- setdiff(which(crit1 & crit2), snp)
  ## Compute R2 for candidate SNPs.
  R2 <- suppressWarnings(cor(as.matrix(gData$markers[, candidateSnps]), gData$markers[, snp]) ^ 2)
  ## Select SNPs based on R2.
  candidateSnpsNames <- names(which(R2[, 1] > minR2))
  return(which(rownames(gData$map) %in% candidateSnpsNames))
}
