#' Plot function for the class \code{GWAS}
#'
#' Creates a plot of an object of S3 class \code{GWAS}. Three types of plots can be made:
#' \itemize{
#' \item{a manhattan plot using the function \code{\link{manhattanPlot}}}
#' \item{a qq plot using the function \code{\link{qqPlot}}}
#' \item{a qtl plot using the function \code{\link{qtlPlot}}}
#' }
#'
#' When making a manhattan plot all markers in the \code{GWAResult} data.frame in the \code{GWAS} object
#' are plotted. Significant SNPs are taken from the \code{signSnp} data.frame and marked in the plot. Also
#' the LOD-threshold is taken from the \code{GWAS} object and plotted. Both \code{signSnp} and \code{thr}
#' can be left empty and will then be ignored in the plot.\cr
#' \code{...} can be used to pass extra arguments to the actual plotting functions. See those fuctions
#' for details.
#'
#' @param x object of class \code{GWAS}.
#' @param ... further arguments to be passed on to the actual plotting functions.
#' @param type string indicating the type of plot to be made. One of "manhattan", "qq" and "qtl".
#' @param trait a string indicating for which trait the results should be plotted. For \code{type}
#' "qtl" all traits are plotted. If \code{x} only contains result for one trait \code{trait} may be
#' \code{NULL}.
#'
#' @seealso \code{\link{manhattanPlot}}, \code{\link{qqPlot}}, \code{\link{qtlPlot}}
#'
#' @export
plot.GWAS <- function(x, ..., type = "manhattan", trait = NULL) {
  type <- match.arg(type, choices = c("manhattan", "qq", "qtl"))
  GWAResult <- x$GWAResult[[1]]
  signSnp <- x$signSnp[[1]]
  if (type != "qtl") {
    if (is.null(trait)) {
      trait <- unique(GWAResult$trait)
      if (length(trait) > 1) {
        stop("Trait not supplied but multiple traits detected in data.")
      }
    } else {
      GWAResult <- GWAResult[GWAResult$trait == trait, ]
      signSnp <- signSnp[signSnp$trait == trait, ]
    }
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
  } else if (type == "qtl") {
    qtlPlot(data = signSnp,
      map = GWAResult[c("chr", "pos")])
  }
}
