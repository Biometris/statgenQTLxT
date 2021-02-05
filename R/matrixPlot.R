#' Create a matrix plot
#'
#' Create a matrix plot of effect sizes and directions per parent for
#' significant snps.
#'
#' @import grDevices graphics
#'
#' @keywords internal
matrixPlot <- function(effectDat,
                       signSnp,
                       map,
                       xLab = "Chromosomes",
                       yLab = "trait",
                       chrBoundaries,
                       founders,
                       ...,
                       output = TRUE) {
  ## Extract central chromosome postions from map.
  ## Differentiate cases to deal with character chromosomes.
  if (is.numeric(map$chr)) {
    chrs <- as.numeric(levels(as.factor(map$chr)))
  } else {
    chrs <- levels(as.factor(map$chr))
  }
  xMarks <- aggregate(x = map$cumPos, by = list(map$chr), FUN = function(x) {
    min(x) + (max(x) - min(x)) / 2
  })[, 2]
  ## Compute chromosome boundaries.
  chrBnd <- cumsum(chrBoundaries[, 2])[-nrow(chrBoundaries)]
  ## Add cumulative position from map to effects.
  parEffData <- merge(effectDat, map[, c("snp", "cumPos")], by = "snp")
  ## Add 5 to cumPos since geom_tile uses center of tile and has width 10.
  parEffData$cumPos <- parEffData$cumPos + 5
  ## Only plotting the effects for significant SNPs. Set all others to NA.
  parEffData[!parEffData$snp %in% signSnp$snp, "effect"] <- NA
  p <- ggplot2::ggplot(data = parEffData,
                       ggplot2::aes_string(x = "cumPos", y = "trait",
                                           fill = "effect")) +
    ggplot2::geom_tile(ggplot2::aes(height = 0.5, width = 10)) +
    ggplot2::scale_fill_gradientn(colors = c("blue", "cyan", "white",
                                             "yellow","red"),
                                  values = scales::rescale(c(min(parEffData$effect, na.rm = TRUE),
                                                             0 - sqrt(.Machine$double.eps),
                                                             0, 0 + sqrt(.Machine$double.eps),
                                            max(parEffData$effect, na.rm = TRUE))),
                                  na.value = "white") +
    ggplot2::scale_x_continuous(breaks = xMarks, labels = chrs,
                                expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::geom_vline(xintercept = chrBnd, color = "grey20",
                        lty = 2, size = 0.3) +
    ggplot2::labs(x = xLab, y = yLab) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        color = "black",
                                                        size = 0.5,
                                                        linetype = "solid"))
  if (output) {
    plot(p)
  }
  invisible(p)
}


