#' Create a manhattan plot
#'
#' Given vectors of marker positions and corresponding LOD-scores plot a
#' LOD-profile. Significant markers can be highlighted with red dots. If there
#' are previously known marker effect also false positives and true negatives
#' will be marked.
#'
#' @param xValues A vector of cumulative marker positions.
#' @param yValues A vector of LOD-scores.
#' @param map A data.frame with at least the columns chr, the number of the
#' chromosome and cumPos, the cumulative position of the snp the cumulative
#' position of the snp starting from the first chromosome.
#' @param fileName A character string, the name of the outputfile that is
#' created. If left empty the plot is written to the screen.
#' @param jpegPlot Should a jpeg file be produced? If \code{FALSE} a pdf file
#' is produced.
#' @param xLab A character string, the x-axis label.
#' @param yLab A character string, the y-axis label.
#' @param plotType The type of plot that is made, either a line plot ("l") or
#' a dot plot ("d" or "p")
#' @param xSig A vector of integers, indicating which components in the vectors
#' xValues and yValues are significant.
#' @param xEffects A vector of integers, indicating which components in the
#' vector xValues correspond to a real (known) effect.
#' @param colPalette A color palette used for plotting.
#' @param chrBoundaries A vector of chromosome boundaries, i.e. x-values on the
#' same scale as xValues.
#' @param yThr A numerical value for the LOD-threshold.
#' @param signPointsThickness A numerical value giving the thickness of the
#' points that are false/true positives/negatives.
#' @param ... Other graphical parameters passed on to the actual plot function.
#'
#' @return A LOD-profile with LOD-scores per snip. Markers declared significant
#' get a red dot, markers with a real effect get a blue dot. If both significant
#' and real effects are given false positives get an orange dot, true negatives
#' a yellow dot and true positives a green dot.
#'
#' @import grDevices graphics
#'
#' @export
manhattanPlot <- function(xValues,
                          yValues,
                          map,
                          fileName = "",
                          jpegPlot = TRUE,
                          xLab = "Chromosomes",
                          yLab = expression(-log[10](p)),
                          plotType = c("l", "p", "d"),
                          xSig = integer(),
                          xEffects = integer(),
                          colPalette = rep(c("royalblue", "maroon"), 50),
                          chrBoundaries = 0,
                          yThr = NULL,
                          signPointsThickness = 0.6,
                          ...) {
  ## Basic argument checks
  if (is.null(xValues) || !is.numeric(xValues)) {
    stop("xValues should be an integer vector")
  }
  if (is.null(yValues) || !is.numeric(yValues)) {
    stop("yValues should be a numerical vector")
  }
  if (fileName != "" && (is.null(fileName) || length(fileName) > 1 ||
                         !is.character(fileName))) {
    stop("fileName cannot be empty")
  }
  if (fileName != "" && (is.null(jpegPlot) || length(jpegPlot) > 1 ||
                         !is.logical(jpegPlot))) {
    stop("jpegPlot should be a single logical")
  }
  if (is.null(xSig) || !is.numeric(xSig)) {
    stop("xSig should be an integer vector")
  }
  if (is.null(xEffects) || !is.numeric(xEffects)) {
    stop("xEffects should be an integer vector")
  }
  if (is.null(chrBoundaries) || !is.numeric(chrBoundaries)) {
    stop("chrBoundaries should be an integer vector")
  }
  if (!is.null(yThr) && (length(yThr) > 1)) {
    stop("yThr should be a single integer")
  }
  if (is.null(signPointsThickness) || length(signPointsThickness) > 1 ||
      !is.numeric(signPointsThickness)) {
    stop("signPointsThickness should be a single number")
  }
  ## Check correspondence xValues and yValues
  if (length(xValues) != length(yValues)) {
    stop("xValues and yValues should be of the same length")
  }
  if (!(plotType %in% c("l", "d", "p")) ) {
    plotType <- "l"
  }
  ## Open file connection
  if (fileName != "") {
    if (jpegPlot) {
      jpeg(fileName, width = 720, height = 480, quality = 100)
    } else {
      pdf(fileName)
    }
  }
  ## Extract central chromosome postions from map
  ## Differentiate cases to deal with character chromosomes
  if (is.numeric(map$chr)) {
    chromosomes <- as.numeric(levels(as.factor(map$chr)))
  } else {
    chromosomes <- levels(as.factor(map$chr))
  }
  xMarks <- aggregate(x = map$cumPos, by = list(map$chr),
                      FUN = function(x) {
                        min(x) + (max(x) - min(x)) / 2
                      })[, 2]
  ## Setup empty plot
  ## ylim has to be set to cope with subsetting based on lod.
  plot(x = xValues, y = yValues, ylim = c(0, max(yValues, na.rm = TRUE)),
       xlab = xLab, ylab = yLab, type = "n", lwd = 0.4,
       xaxt = 'n', ...)
  axis(side = 1, at = xMarks, labels = chromosomes, cex.axis = 0.8)
  ## If chromosome boundaries are known add lines/ points per chromosome
  if (sum(chrBoundaries) != 0) {
    for (chromosome in chromosomes) {
      if (plotType == "l") {
        lines(x = xValues[map$chr == chromosome],
              y = yValues[map$chr == chromosome],
              lwd = 0.4, col = colPalette[which(chromosomes == chromosome)])
      } else {
        points(x = xValues[map$chr == chromosome],
               y = yValues[map$chr == chromosome], pch = 20, lwd = 0.4,
               col = colPalette[which(chromosomes == chromosome)])
      }
    }
  } else {
    ## If chromosome boundaries are unknown add lines/ points
    if (plotType == "l") {
      lines(x = xValues, y = yValues, xlab = xLab, ylab = yLab, lwd = 0.4,
            col = "royalblue")
    } else {
      points(x = xValues, y = yValues, xlab = xLab, ylab = yLab, pch = 20,
             lwd = 0.4, col = "royalblue")
    }
  }
  ## Add red dots for significant markers
  if (sum(xSig) != 0 && sum(xEffects) == 0) {
    points(x = xValues[xSig], y = yValues[xSig],
           pch = 20, col = "red", lwd = signPointsThickness)
  } else if (sum(xSig) == 0 && sum(xEffects) != 0) {
    ## Add blue dots for known effects
    points(x = xValues[xEffects], y = yValues[xEffects],
           pch = 20,col = "blue", lwd = signPointsThickness)
  } else if (sum(xSig) != 0 && sum(xEffects) != 0) {
    ## Add orange/yellow/green dots for false positives/ true negatives and
    ## true positives
    falsePos <- setdiff(xSig, xEffects)
    trueNeg <- setdiff(xEffects, xSig)
    truePos <- intersect(xSig, xEffects)
    points(x = xValues[falsePos], y = yValues[falsePos],
           pch = 20, col = "orange", lwd = signPointsThickness)
    points(x = xValues[trueNeg], y = yValues[trueNeg],
           pch = 20, col = "yellow", lwd = signPointsThickness)
    points(x = xValues[truePos], y = yValues[truePos],
           pch = 20, col = "green", lwd = signPointsThickness)
  }
  ## Add a horizontal line for the threshold
  if (!is.null(yThr)) {
    abline(h = yThr, lty = 2, lwd = 0.5)
  }
  ## Close file connection
  if (fileName != "") {
    dev.off()
  }
}



