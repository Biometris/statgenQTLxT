#' qq-plot of observed versus expected LOD-scores
#'
#' Given a vector of pvalues, generate a qq-plot of observed LOD-scores versus
#' expected LOD-scores. Code taken from Segura et al. and adapted
#'
#' @inheritParams manhattanPlot
#'
#' @param pValues A numeric vector of pValues. Missings are ignored when
#' plotting.
#'
#' @references Segura et al. (2012) An efficient multi-locus mixed model
#' approach for genome-wide association studies in structured populations
#'
#' @import grDevices graphics
#'
#' @export
qqPlot <- function(pValues,
                   fileName = "",
                   jpegPlot = TRUE,
                   ...,
                   output = TRUE) {
  dotArgs <- list(...)
  if (is.null(pValues) || !is.numeric(pValues) || any(pValues < 0) ||
      any(pValues > 1)) {
    stop("pValues should be an numeric vector with values between 0 and 1")
  }
  if (fileName != "" && (is.null(fileName) || length(fileName) > 1 ||
                         !is.character(fileName))) {
    stop("fileName cannot be empty")
  }
  if (fileName != "" && (is.null(jpegPlot) || length(jpegPlot) > 1 ||
                         !is.logical(jpegPlot))) {
    stop("jpegPlot should be a single logical")
  }
  if (fileName != "") {
    if (jpegPlot) {
      jpeg(fileName, width = 720, height = 480, quality = 100)
    } else {
      pdf(fileName)
    }
  }
  pValues <- na.omit(pValues)
  expected <- -log10(ppoints(n = length(pValues)))
  observed <- -log10(sort(pValues))
  plotDat <- data.frame(expected, observed)
  ## Construct title.
  if (!is.null(dotArgs$title)) {
    plotTitle <- dotArgs$title
  } else {
    plotTitle <- "QQ-plot"
  }
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes_string(x = "expected", y = "observed")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = expression(Expected~~-log[10](p)),
                  y = expression(Observed~~-log[10](p))) +
    ggplot2::ggtitle(plotTitle)
  if (output) {
    plot(p)
  }
  if (fileName != "") {
    dev.off()
  }
  invisible(p)
}
