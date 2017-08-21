#' Create a GWAS object
#'
#' \code{createGWAS} creates an object of S3 class containing the results of a GWAS analysis.
#' \code{GWAResult} and \code{signSnp} are both optional, however at least one of those should be provided
#' as input.
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
#'
#' @return an object of class GWAS, a list of the input items.
#'
#' @return A list containing two data.frames:
#' \code{GWAResult}, the full results for all markers with the following columns:
#' \itemize{
#' \item{trait: trait name.}
#' \item{snp: marker name.}
#' \item{chr: chromosome on which the marker lies.}
#' \item{pos: position of the marker on the chromosome.}
#' \item{pValue: p-value of the GLS-test.}
#' \item{effect: effect size.}
#' \item{effectSe: standard error of the effect size.}
#' \item{LOD: LOD-score.}
#' \item{RLR2: likelihood-ratio based \eqn{R^2} as described by Sun et al.}
#' \item{allFreq: allele frequency.}
#' }
#' \code{SignSnp}, results for significant SNPs including SNPs close to significant SNPs if
#' \code{sizeInclRegion} > 0. \code{SignSnp} contains the following columns:
#' \itemize{
#' \item{trait: trait name.}
#' \item{snp: marker name.}
#' \item{chr: chromosome on which the marker lies.}
#' \item{pos: position of the marker on the chromosome.}
#' \item{pValue: p-value of the GLS-test.}
#' \item{LOD: LOD-score.}
#' \item{snpStatus: status of the SNP, i.e. "significant SNP" or "within ... of sign. SNP".}
#' \item{allFreq: allele frequency.}
#' \item{effect: effect size.}
#' \item{RLR2: likelihood-ratio based \eqn{R^2} as described by Sun et al.}
#' \item{propSnpVar: proportion of the variance explained by the SNP.}
#' }
#'
#'

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
  ## Check signSnps
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

is.GWAS <- function(x) {
  inherits(x, "GWAS")
}

summary.GWAS <- function(object) {
  GWAResult <- object$GWAResult[[1]]
  signSnp <- object$signSnp[[1]]
  GWASInfo <- object$GWASInfo
  ## Print traits
  cat("Traits analysed: ", paste(unique(GWAResult$trait), collapse = ", "), "\n\n")
  ## Print SNP numbers.
  cat("Data are available for", nrow(GWAResult), "SNPs.\n")
  if (!is.null(GWASInfo$MAF)) {
    cat(sum(is.na(GWAResult$pValue)), "of them were not analyzed because their minor allele frequency",
      "is below", GWASInfo$MAF, "\n\n")
  }
  if (as.numeric(GWASInfo$GLSMethod) == 1) {
    ## Print mixed model info.
    cat("Mixed model with only polygenic effects, and no marker effects:\n")
    cat("Genetic variance: ", GWASInfo$varcomp[1], "\n")
    cat("Residual variance: ", GWASInfo$varcomp[2], "\n\n")
  }
  if (as.numeric(GWASInfo$thrType) %in% 1:3) {
    ## Print significant SNP info.
    cat("LOD-threshold: ", object$thr, "\n")
    if (!is.null(signSnp)) {
      cat("Number of significant SNPs =" , nrow(signSnp[signSnp$snpStatus == "significant snp", ]), "\n")
      cat("Smallest p-value among the significant SNPs:",
        min(signSnp[signSnp$snpStatus == "significant snp", "pValue"]), "\n")
      cat("Largest  p-value among the significant SNPs:",
        max(signSnp[signSnp$snpStatus == "significant snp", "pValue"]),
        "(LOD-score:", min(signSnp[signSnp$snpStatus == "significant snp", "LOD"]), ")\n")
    } else {
      cat("No significant SNPs found.","\n")
    }
    if (GWASInfo$genomicControl) {
      ## Print genomic control.
      cat("\nGenomic control correction was applied\n")
      cat("Genomic control inflation-factor = ", GWASInfo$inflationFactor, "\n\n")
    } else {
      cat("\n","No Genomic control correction was applied", "\n")
    }
  }
}

plot.GWAS <- function(x, y, ..., type = "manhattan") {
  GWAResult <- x$GWAResult[[1]]
  signSnp <- x$signSnp[[1]]
  if (!is.character(type) || length(type) > 1 || !type %in% c("manhattan", "qq")) {
    stop("type should be one of manhattan and qq.")
  }
  if (type == "manhattan") {
    ## Compute chromosome boundaries.
    chrBnd <- aggregate(GWAResult$pos, by = list(GWAResult$chr), FUN = max)
    ## Compute cumulative positions.
    addPos <- data.frame(chr = chrBnd[, 1], add = c(0, cumsum(chrBnd[,2]))[1:nrow(chrBnd)])
    map <- merge(data.frame(chr = GWAResult$chr, pos = GWAResult$pos), addPos, by = "chr")
    map <- data.frame(snp = GWAResult$snp, chr = map$chr, cumPos = map$pos + map$add)
    ## Extract numbers of significant SNPs.
    signSnpNr <- which(map$snp %in% signSnp$snp[signSnp$snpStatus == "significant snp"])
    ## Create manhattan plot.
    manhattanPlot(xValues = map$cumPos, yValues = GWAResult$LOD, map = map,
      plotType = "p", xSig = signSnpNr, chrBoundaries = chrBnd[ , 2], yThr = x$thr, ...)
  } else if (type == "qq") {
    ## Create qq-plot
    qqPlot(pValues = na.omit(GWAResult$pValue), ...)
  }
}


