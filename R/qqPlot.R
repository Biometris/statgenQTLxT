#' qq-plot of observed versus expected LOD-scores
#'
#' Given a vector of pvalues, generate a qq-plot of observed LOD-scores versus
#' expected LOD-scores. Code taken from Segura et al. and adapted
#'
#' @inheritParams manhattanPlot
#'
#' @param pValues a numeric vector of pValues. Missings are ignored when
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
                   ...) {
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
  pMax <- ceiling(max(observed))
  if (!is.null(dotArgs$main)) {
    main <- dotArgs$main
  } else {
    main <- "QQ-plot"
  }
  if (!is.null(dotArgs$col)) {
    color <- dotArgs$col
  } else {
    color <- 1
  }
  if (!is.null(dotArgs$cex)) {
    cex <- dotArgs$cex
  } else {
    cex <- 0.9
  }
  if (!is.null(dotArgs$pch)) {
    pch <- dotArgs$pch
  } else {
    pch <- 20
  }
  if (!is.null(dotArgs$xlab)) {
    xlab <- dotArgs$xlab
  } else {
    xlab <- expression(Expected~~-log[10](p))
  }
  if (!is.null(dotArgs$ylab)) {
    ylab <- dotArgs$ylab
  } else {
    ylab <- expression(Observed~~-log[10](p))
  }
  do.call(plot, c(list(x = expected, y = observed, type = 'b', pch = pch,
                       cex = cex, col = color, xlab = xlab, ylab = ylab,
                       xlim = c(0, max(expected) + 1), ylim = c(0, pMax),
                       main = main),
                  dotArgs[!names(dotArgs) %in% c("pch", "cex", "col", "main",
                                                 "xlab", "ylab")]))
  abline(a = 0, b = 1, col = "blue")
  if (fileName != "") {
    dev.off()
  }
}
