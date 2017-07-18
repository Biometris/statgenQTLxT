#' qq-plot of observed versus expected LOD-scores
#'
#' Given a vector of pvalues, generate a qq-plot of observed LOD-scores versus expected LOD-scores.
#' Code taken from Segura et al. and adapted
#'
#' @inheritParams lodPlot
#'
#' @param pValues a numeric vector of pValues. Missings are ignored when plotting.
#' @param title main title for the plot
#'
#' @references Segura et al. (2012) An efficient multi-locus mixed model approach for genome-wide
#' association studies in structured populations
#'
#' @export

## TO DO: example

qqPlot <- function(pValues,
  title = "QQ-plot",
  fileName = "",
  jpegPlot = TRUE) {
  pValues <- na.omit(pValues)
  expected <- -log10(ppoints(n = length(pValues)))
  observed <- -log10(sort(pValues))
  pMax <- ceiling(max(observed))

  if (fileName != "") {
    if (jpegPlot) {jpeg(fileName, width = 720, height = 480, quality = 100)} else {pdf(fileName)}
  }

  plot(x = expected, y = observed, type = 'b', pch = 20, cex = 0.9, col = 1,
    xlab = expression(Expected~~-log[10](p)),
    ylab = expression(Observed~~-log[10](p)),
    xlim=c(0, max(expected) + 1),
    ylim=c(0, pMax),
    main = title)
  abline(a = 0, b = 1, col="blue")

  if (fileName != "") dev.off()
}
