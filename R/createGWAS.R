#' S3 Class GWAS
#'
#' \code{createGWAS} creates an object of S3 class GWAS containing the results of a GWAS analysis.
#' \code{GWAResult} and \code{signSnp} are both optional, however at least one of those should be provided
#' as input.\cr
#' \code{summary} and \code{plot} functions are available.\cr\cr
#' \code{is.gData} tests if an \code{R} object is of class \code{gData}.
#'
#' @param GWAResult optional data.frame or list of data.frames containing the overall analysis results.
#' Should a least contain columns \code{trait}, the evaluated trait, \code{snp}, the name of the SNP,
#' \code{chr}, the chromosome number, \code{pos}, the position of the SNP on the chromosome,
#' \code{pValue}, the p-values from the analysis and \code{LOD} the LOD-score.
#' @param signSnp optional data.frame or list of data.frames containing information on the significant
#' SNPs and optionally the SNPs close to the significant SNPs. Should at least contain columns
#' \code{trait}, the evaluated trait, \code{snp}, the name of the SNP, \code{pValue}, the p-values
#' from the analysis and \code{LOD} the LOD-score.
#' @param kin optional kinship matrix or list of chromosome specific kinship matrices.
#' @param thr optional numeric value, the threshold used in performing the GWAS analysis.
#' @param GWASInfo list containing extra information concering the GWAS analysis.
#' @param x \code{R} object
#'
#' @return \code{createGWAS} returns an object of class GWAS, a list of the input items.\cr\cr
#' \code{is.gData} returns \code{TRUE} or \code{FALSE} depending on whether its argument is a \code{gData}
#' object.
#'
#' @seealso \code{\link{summary.GWAS}}, \code{\link{plot.GWAS}}
#'
#' @name GWAS
NULL

#' @rdname GWAS
#' @export
createGWAS <- function(GWAResult = NULL,
  signSnp = NULL,
  kin = NULL,
  thr = NULL,
  GWASInfo = NULL) {
  ## Check that at least one input is provided.
  if (is.null(GWAResult) && is.null(signSnp))
    stop("At least one of GWAResult and signSnp should be provided.\n")
  ## Check GWAResults
  if (!is.null(GWAResult)) {
    if (!is.data.frame(GWAResult) &&
        !(is.list(GWAResult) && all(sapply(GWAResult, FUN = is.data.frame))))
      stop("GWAResult should be a data.frame or a list data.frames.\n")
    if (is.data.frame(GWAResult)) {
      ## If not a list already put data.frame in a list.
      GWAResult <- list(GWAResult)
    }
    if (!all(sapply(GWAResult, FUN = function(x) {
      all(c("trait", "snp", "chr", "pos", "pValue", "LOD") %in% colnames(x))})))
      stop("GWAResult should contain columns trait, snp, chr, pos, pValue and LOD.\n")
  }
  ## Check signSnps
  if (!is.null(signSnp)) {
    if (!is.data.frame(signSnp) &&
        !(is.list(signSnp) && all(sapply(signSnp, FUN = is.data.frame))))
      stop("signSnp should be a data.frame or a list of data.frames.\n")
    if (is.data.frame(signSnp)) {
      ## If not a list already put data.frame in a list.
      signSnp <- list(signSnp)
    }
    if (!all(sapply(signSnp, FUN = function(x) {
      all(c("trait", "snp", "snpStatus", "pValue", "LOD") %in% colnames(x))})))
      stop("signSnp should contain columns trait, snp, snpStatus, pValue and LOD.\n")
  }
  ## Check kin
  if (!is.null(kin)) {
    if (!is.matrix(kin) &&
        !(is.list(kin) && all(sapply(kin, FUN = is.matrix))))
      stop("kin should be a matrix or a list of matrices.\n")
  }
  ## Create GWAS object.
  GWAS <- structure(list(GWAResult = GWAResult,
    signSnp = signSnp,
    kinship = kin,
    thr = thr,
    GWASInfo = GWASInfo),
    class = "GWAS")
  return(GWAS)
}

#' @rdname GWAS
#' @export
is.GWAS <- function(x) {
  inherits(x, "GWAS")
}








