matrixPlot <- function(data,
                       signSnp) {
  MPQoutput <- data
  nchr = nlevels(factor(MPQoutput$chr))
  cmin = tapply(MPQoutput$pos, MPQoutput$chr, min)
  cmax = tapply(MPQoutput$pos, MPQoutput$chr, max)
  range = (cmax - cmin) + 8
  sumrange = cumsum(range)
  start_chr = c(0, sumrange[1:(nchr - 1)])
  ver_line = sumrange[1:(nchr - 1)]
  xtics = (start_chr + sumrange) * 0.5
  x = start_chr[MPQoutput$chr] + (MPQoutput$pos - cmin[MPQoutput$chr]) + 4
  parEffData <- MPQoutput %>% mutate(x = x) %>%
    gather('parent','effect', grep('eff.', names(MPQoutput)))
  QTL.parEffData <- filter(parEffData) %>% filter(name %in% signSnp$name)
  parEffData$parent <- factor(parEffData$parent,
                              levels = unique(parEffData$parent))
  p <- ggplot2::ggplot(parEffData, ggplot2::aes_string("x", "parent")) +
    ggplot2::geom_tile(ggplot2::aes_string(fill = "effect", height = 1,
                                           width = 10)) +
    ggplot2::scale_fill_gradient2(low = "dodgerblue", mid = 'white',
                                  high = "red", midpoint = 0, space = "Lab",
                                  guide = "colourbar") +
    ggplot2::geom_vline(xintercept = ver_line, color = "grey20",
                        lty = 2, size = 0.3) +
    ggplot2::labs(x = "Chromosome", y = "Parent") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        color = "black",
                                                        size = 0.5,
                                                        linetype = "solid"))
  return(p)
}
