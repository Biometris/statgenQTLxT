#' recoding markers and imputation
#'
#' recoding markers and imputation
#'
#' @param gData an object of class \code{gData} containing at least \code{markers}.
#' @param nMiss a numical value between 0 and 1. SNPs with a fraction of missings higher then
#' \code{nMiss} will be removed.
#' @param keep a vector of SNPs. These SNPs will never be removed in the whole process.

recodeMarkers <- function(gData,
                          nMiss = NULL,
                          keep = NULL) {
  markersOrig <- gData$markers


  ## TEST ONLY
  markersOrig <- matrix(c(
    "AA",   "AA",   "AA",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
    "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",  NA,
    "AA",   "AA",   "AB",   "BB",   "AB",   "AA",   "AA",   "BB",  NA,
    "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",  NA,
    "AA",   "AA",   "BB",   "AB",   "AA",   "BB",   "BB",   "BB",  "AB",
    "AA",   "AA",   "BB",   "BB",   "AA",   NA,     "BB",   "AA",  NA,
    "AB",   "AA",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",  NA,
    "AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA",
    "AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA",
    "AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",  NA),
    ncol=9,byrow=TRUE)
  colnames(markersOrig) <- paste("SNP",1:9,sep="")
  rownames(markersOrig) <- paste("ID",1:10+100,sep="")

  snpKeep <- colnames(markersOrig) %in% keep
  ## Remove markers with too many missings
  if (!is.null(nMiss)) {
    snpMiss <- colMeans(is.na(markersOrig)) <= nMiss
    markersNew <- markersOrig[, snpMiss | snpKeep]
  }
}
