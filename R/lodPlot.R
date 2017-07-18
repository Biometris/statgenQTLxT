#' Plot a LOD-profile
#'
#' Given vectors of marker positions and corresponding LOD-scores plot a LOD-profile.
#' Significant markers can be highlighted with red dots. If there are previously known
#' marker effect also false positives and true negatives will be marked.
#'
#' @param xValue vector of cumulative marker positions
#' @param yValues vector of LOD-scores
#' @param map A dataframe with four columns; chromosome, the number of the chromosome,
#' position, the position of the snp on the chromosome, snp.name, the name of the snp and
#' cum.position, the cumulative position of the snp starting from the first chromosome.
#' @param fileName name of the outputfile that is created. If left empty the plot is written
#' to the screen.
#' @param jpegPlot should a jpeg file be produced? If \code{FALSE} a pdf file is produced.
#' @param xLab x-axis label
#' @param yLab y-axis label
#' @param title title of the plot
#' @param plotType lines ("l") or dots ("d" or "p")
#' @param xSig vector of integers, indicating which components in the vectors xValues and
#' yValues are significant.
#' @param xEffects vector of integers, indicating which components in the vector xValues
#' correspond to a real (known) effect.
#' @param colPalette color palette used for plotting.
#' @param chrBoundaries vector of chromosome boundaries, i.e. x-values on the same scale as xValues.
#' @param yThr LOD-threshold.
#' @param signPointsThickness thickness of the points that are false/true positives/negatives.
#'
#' @return a LOD-profile with LOD-scores per snip. Markers declared significant get a red dot,
#' markers with a real effect get a blue dot. If both significant and real effects are given
#' false positives get an orange dot, true negatives a yellow dot and true positives a green dot.
#'
#' @export

## TO DO: check input
## TO DO: example

lodPlot <- function(xValues,
  yValues,
  map,
  fileName = "",
  jpegPlot = TRUE,
  xLab = "Chromosomes",
  yLab = expression(-log[10](p)),
  title = "",
  plotType = "l",
  xSig = integer(),
  xEffects = integer(),
  colPalette = rep(c("royalblue", "maroon"), 50)[1:length(levels(factor(map$chromosome)))],
  chrBoundaries = 0,
  yThr = 0,
  signPointsThickness = 0.6) {

  if (!(plotType %in% c("l","d","p")) ) {plotType <- "l"}
  if (plotType == "p") {plotType <- "d"}

  if (fileName != "") {
    if (jpegPlot) {jpeg(fileName, width = 720, height = 480, quality = 100)} else {pdf(fileName)}
  }

  chromosomes <- as.numeric(levels(factor(map$chromosome)))
  xMarks <- aggregate(x = map$cum.position, by = list(map$chromosome),
    FUN = function(x) {min(x) + (max(x) - min(x)) / 2})[, 2]
  plot(x = xValues, y = yValues, xlab = xLab, ylab = yLab, type = "n", lwd = 0.4,
    main = title, xaxt = 'n')
  axis(side = 1, at = xMarks, labels = chromosomes)

  if (sum(chrBoundaries) != 0) {
    for (chromosome in chromosomes) {
      if (plotType == "l") {
        lines(x = xValues[map$chromosome == chromosome],
          y = yValues[map$chromosome == chromosome],
          type = "l", lwd = 0.4, col = colPalette[which(chromosomes == chromosome)])
      } else {
        points(x = xValues[map$chromosome == chromosome],
          y = yValues[map$chromosome == chromosome], pch=20, lwd=0.4,
          col = colPalette[which(chromosomes==chromosome)])
      }
    }
  } else {
    if (plotType == "l") {
      lines(x = xValues, y = yValues, xlab = xLab, ylab = yLab, lwd = 0.4, col = "royalblue")
    } else {
      points(x = xValues, y = yValues, xlab = xLab, ylab = yLab, pch = 20, lwd = 0.4, col = "royalblue")
    }
  }
  if (sum(xSig) != 0 && sum(xEffects) == 0) {
    points(x = xValues[xSig], y = yValues[xSig],
      pch = 20, col = "red", lwd = signPointsThickness)
  } else if (sum(xSig) == 0 && sum(xEffects) != 0) {
    points(x = xValues[xEffects], y = yValues[xEffects],
      pch = 20,col = "blue", lwd = signPointsThickness)
  } else if (sum(xSig) != 0 && sum(xEffects) != 0) {
    falseNeg <- setdiff(xEffects, xSig)
    falsePos <- setdiff(xSig, xEffects)
    truePos <- intersect(xSig, xEffects)
    points(x = xValues[falsePos], y = yValues[falsePos],
      pch = 20, col = "orange", lwd = signPointsThickness)
    points(x = xValues[falseNeg], y = yValues[falseNeg],
      pch = 20, col = "yellow", lwd = signPointsThickness)
    points(x = xValues[truePos], y = yValues[truePos],
      pch = 20, col = "green", lwd = signPointsThickness)
  }
  if (yThr != 0) abline(h = yThr, lty = 2, lwd = 0.5)
  if (fileName != "") dev.off()
}



