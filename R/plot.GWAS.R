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
#' are plotted. Significant SNPs are taken from the \code{signSnp} data.frame in the \code{GWAS} object and
#' marked in the plot. Also the LOD-threshold is taken from the \code{GWAS} object and plotted as
#' a horizontal line. Both \code{signSnp} and \code{thr} can be left empty and will then be
#' ignored in the plot.\cr
#' \code{...} can be used to pass extra arguments to the actual plotting functions. See those fuctions
#' for details.
#'
#' @param x object of class \code{GWAS}.
#' @param ... further arguments to be passed on to the actual plotting functions.
#' @param type string indicating the type of plot to be made. One of "manhattan", "qq" and "qtl".
#' @param environment a string or numeric index indicating for which environment the
#' plot should be made. If \code{x} only contains results for one trait \code{environment} may be
#' \code{NULL}.
#' @param trait a string indicating for which trait the results should be plotted. For \code{type}
#' "qtl" all traits are plotted. If \code{x} only contains results for one trait \code{trait} may be
#' \code{NULL}.
#'
#' @seealso \code{\link{manhattanPlot}}, \code{\link{qqPlot}}, \code{\link{qtlPlot}}
#'
#' @export

plot.GWAS <- function(x, ..., type = "manhattan", environment = NULL, trait = NULL) {
  type <- match.arg(type, choices = c("manhattan", "qq", "qtl"))
  ## Checks.
  if (!is.null(environment) && !is.character(environment) && !is.numeric(environment))
    stop("environment should be a character or numeric value.\n")
  if ((is.character(environment) && !environment %in% names(x$GWAResult)) ||
      (is.numeric(environment) && !environment %in% 1:length(x$GWAResult)))
    stop("environment should be in x.\n")
  ## Convert character input to numeric.
  if (is.character(environment)) environment <- which(names(x$GWAResult) == environment)
  ## If NULL then summary of all environment.
  if (is.null(environment)) {
    if (length(x$GWAResult) != 1) {
      stop("Environment not supplied but multiple environments detected in data.")
    } else {
      environment <- 1
    }
  }
  GWAResult <- x$GWAResult[[environment]]
  signSnp <- x$signSnp[[environment]]
  if (type != "qtl") {
    if (is.null(trait)) {
      trait <- unique(GWAResult$trait)
      if (length(trait) > 1) {
        if (substr(as.character(x$GWASInfo$call)[1], 1, 9) == "runSingle") {
          stop("Trait not supplied but multiple traits detected in data.")
        } else {
          ## For multi trait GWAS p-values are the same for all traits.
          trait <- trait[1]
        }
      }
    } else {
      GWAResult <- GWAResult[GWAResult$trait == trait, ]
      signSnp <- signSnp[signSnp$trait == trait, ]
    }
  }
  if (type == "manhattan") {
    ## Compute chromosome boundaries.
    GWAResult <- GWAResult[!is.na(GWAResult$pos), ]
    chrBnd <- aggregate(GWAResult$pos, by = list(GWAResult$chr), FUN = max)
    ## Compute cumulative positions.
    addPos <- data.frame(chr = chrBnd[, 1], add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)])
    map <- merge(data.frame(chr = GWAResult$chr, pos = GWAResult$pos), addPos, by = "chr")
    map <- data.frame(snp = GWAResult$snp, chr = map$chr, cumPos = map$pos + map$add)
    ## Extract numbers of significant SNPs.
    signSnpNr <- which(map$snp %in% signSnp$snp[as.numeric(signSnp$snpStatus) == 1])
    ## Create manhattan plot.
    manhattanPlot(xValues = map$cumPos, yValues = GWAResult$LOD, map = map,
      plotType = "p", xSig = signSnpNr, chrBoundaries = chrBnd[, 2], yThr = x$thr, ...)
  } else if (type == "qq") {
    ## Create qq-plot
    qqPlot(pValues = na.omit(GWAResult$pValue), ...)
  } else if (type == "qtl") {
    qtlPlot(data = signSnp,
      map = GWAResult[!is.na(GWAResult$pos), c("chr", "pos")])
  }
}
