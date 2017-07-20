#' Draw a qtl plot
#'
#' Description to be added.
#'
#' @param data QTL data to be plotted
#' @param chromosome a character indicating the data column containing the chromosome number.
#' @param trait a character indicating the data column containing the trait name.
#' @param snpEffect a character indicating the data column containing the snp effect.
#' @param snpPosition a character indicating the data column containing the position of the snp on
#' the chromosome.
#' @param map a dataframe with at least the columns chromosome, the number of the chromosome
#' and position, the position of the snp on the chromosome. These are used for calculating the
#' fysical limits of the chromosome.
#' @param normalize should the snpEffect be normalized?
#' @param sortData should the data be sorted before plotting? Either false if no sorting should be done
#' or a character indicating the data column on which the data should be sorted.
#' @param yLab y-axis label.
#' @param exportPptx should the plot be exported to a .pptx file.
#' @param pptxName name of the .pptx file to which the plot is exported. Ignored if exportPptx =
#' \code{FALSE}
#'
#' @return a plot.
#'
#' @references Millet et al. (2016) Genome-wide analysis of yield in Europe: Allelic effects vary
#' with drought and heat scenarios
#'
#' @export


## TO DO: example
## # Example 1 : different phenotypic trait
##snpList <- read.table("./example_data_drops/ListeQTL_LOD45_DROPS_Response.csv",sep=",",header=TRUE)
# Example 2: different environment
##snpList2 <- read.table("./example_data_drops/Fixed_modelMEML_33+15QTL_GY_perenvt_effB73.csv",sep=",",header=TRUE)
## TO DO: bin positions
## TO DO: description + output
## Dataframe format should be output of the "listQTL" function

qtlPlot <- function(data,
  chromosome = "chromosome",
  trait = "trait",
  snpEffect = "effect",
  snpPosition = "position",
  map,
  normalize = FALSE,
  sortData = FALSE,
  yLab = "Traits",
  exportPptx = FALSE,
  pptxName = "") {

  ## Basic argument checks
  if (is.null(data) || !is.data.frame(data)) stop("data should be a dataframe")
  if (is.null(chromosome) || length(chromosome) > 1 || !is.character(chromosome))
    stop("chromosome should be a single character")
  if (is.null(trait) || length(trait) > 1 || !is.character(trait))
    stop("trait should be a single character")
  if (is.null(snpEffect) || length(snpEffect) > 1 || !is.character(snpEffect))
    stop("snpEffect should be a single character")
  if (is.null(snpPosition) || length(snpPosition) > 1 || !is.character(snpPosition))
    stop("snpPosition should be a single character")
  if (is.null(map) || !is.data.frame(map)) stop("map should be a dataframe")
  if (is.null(normalize) || length(normalize) > 1 || !is.logical(normalize))
    stop("normalize should be a single logical")
  if (is.null(sortData) || (is.logical(sortData) && sortData) ||
      (is.character(sortData) && length(sortData) > 1))
    stop("sortData should be either FALSE or a single character")
  if (is.null(exportPptx) || length(exportPptx) > 1 || !is.logical(exportPptx))
    stop("exportPptx should be a single logical")
  if (exportPptx && (is.null(pptxName) || length(pptxName) > 1 || !is.character(pptxName)))
    stop("pptxName cannot be empty")

  ## Check that all necessary columns are in the data
  requiredColumns <- c(chromosome, trait, snpEffect, snpPosition)
  if (is.character(sortData)) requiredColumns <- c(requiredColumns, sortData)
  requiredCheck <- requiredColumns %in% colnames(data)
  if (!all(requiredCheck))
    stop("data lacks the following columns: ",
      paste0(requiredColumns[!requiredCheck], collapse = ", "), ".\n\n")

  ## Check that all necessary columns are in the map file
  requiredColumnsMap <- c("chromosome", "position")
  requiredCheckMap <- requiredColumnsMap %in% colnames(map)
  if (!all(requiredCheckMap))
    stop("data lacks the following columns: ",
      paste0(requiredColumnsMap[!requiredCheckMap], collapse = ", "), ".\n\n")

  ## Center and reduce the allelic effect (because of the different units)
  if (normalize) {
    data$eff <- sapply(1:nrow(data), function(x)
      (data[x, snpEffect] - mean(data[data[trait] == as.character(data[x, trait]), snpEffect], na.rm=TRUE)) /
        sd(data[data[trait] == as.character(data[x, trait]), snpEffect], na.rm = TRUE) )
  } else data$eff <- data[snpEffect]

  if (is.character(sortData)) {
    data$sort <- data[[sortData]]
  } else {
    data$sort <- 1
  }

  ## Add the physical limits of the chromosomes, calculated from the map file
  ## This ensures plotting of all chromosomes
  limitsLow <- aggregate(map$position, by = list(map$chromosome), FUN = min)
  limitsHigh <- aggregate(map$position, by = list(map$chromosome), FUN = max)
  ## empty dataframe with 2 lines per chromosomes and as many columns as the QTL dataframe
  limits <- data.frame(matrix(ncol = ncol(data), nrow = 2 * nrow(limitsLow)))
  names(limits) <- names(data)
  ## Trait and sort have to be filled. Value is not important
  limits[trait] <- data[1, trait]
  limits$sort <- data[1, "sort"]
  ## Set eff to suppress warnings in printing. Setting it to -Inf ensures nothing is plotted
  limits$eff <- -Inf
  limits[, c(chromosome, snpPosition)] <- rbind(limitsLow, limitsHigh)
  data <- rbind(data, limits)

  ## Select and rename relevant columns for plotting
  plotData <- dplyr::select(data, trait = trait, chromosome = chromosome,
    snpEffect = snpEffect, snpPosition = snpPosition, sort, eff)

  ## Add a column with the allelic effect direction (for points color)
  plotData$color <- ifelse(plotData$eff > 0, "pos", "neg")
  plotData <- droplevels(plotData)

  ### Loading the bin position to add vertical lines to the figure
  ### WARNINGS : particular case of the maize genome !
  ### other species do not have bins
  # dat.vline<-read.table("/home/millet/Documents/PhD-1A/GWAS/A_50K+600K_GyGnbGsAnth/PosBinRefGenV2.csv",sep=",",he=T)
  # names(dat.vline)[1] <- "Chromosome"

  ## Create theme for plot
  emptyThemeQtlPlot <-
    ggplot2::theme(plot.background = ggplot2::element_blank(),
      panel.grid.major= ggplot2::element_line(color = "white"),
      panel.grid.minor= ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 20, face = "bold", vjust = 2) ,
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "gray"),
      axis.text.x = ggplot2::element_text(size = 0),
      axis.text.y = ggplot2::element_text(size = 13 , color = "black"),
      axis.title.x = ggplot2::element_text(size = 20, color = "navyblue"),
      axis.title.y = ggplot2::element_text(size = 20, color = "navyblue"),
      strip.background = ggplot2::element_rect(fill = "gray40"),
      strip.text = ggplot2::element_text(),
      strip.text.x = ggplot2::element_text(size = 14),
      strip.text.y = ggplot2::element_text(size = 0))


  ## Create the plot object
  qtlPlot <-
    ggplot2::ggplot(data = plotData,
      ggplot2::aes(x = snpPosition,
        ## Y data is sorted in reverse order because of the way ggplot plots
        y = reorder(trait, -sort),
        ## Point size proportional to (absolute value of) allelic effect
        size = abs(eff),
        # Point color depends on the effect direction
        color = factor(color)))  +
    ## use custom made theme
    emptyThemeQtlPlot +
    ggplot2::ylab(yLab) +
    ggplot2::xlab("Chromosomes")  +
    # geom_vline(aes(xintercept = position),              # add vertical lines at the bin position
    #   data = dat.vline,                         # the bin positions are contained in the 'dat.vline' dataframe
    #   linetype=1,                               # type of line (1 = full)
    #   color = "white")   +                      # color of line
    ## Add the points with a slight transparency in case of overlap
    ggplot2:: geom_point(alpha = I(0.7)) +
    ## Split of the plot according to the chromosomes on the x axis
    ggplot2::facet_grid(". ~ chromosome",
      ## Do not resize the x axise (otherwise every chromosome has the same size)
      scales = "free",
      ## Do not add extra space between two facets
      space = "free",
      ## Place the chromosome labels at the bottom
      switch = "both") +
    ## Ascribe a color to the allelic direction (column 'color')
    ggplot2::scale_color_manual("color",
      labels = c("neg", "pos"),
      values = c("darkblue", "green4"))

  ## Plot the plot object on screen
  plot(qtlPlot)

  if (exportPptx) {
    ## Save figure in .pptx
    if (requireNamespace("ReporteRs", quietly = TRUE)) {
      ## Create empty .pptx file
      pptOut <- ReporteRs::pptx()
      ## Add new slide (always necessary)
      pptOut <- ReporteRs::addSlide(doc = pptOut, slide.layout = "Title and Content")
      ## Add plot to the document
      pptOut <- ReporteRs::addPlot(doc = pptOut,
        ## "print" the graph
        fun = print,
        x = qtlPlot,
        ## Plot position on the slide in inches
        offx = 1, offy = 1,
        ## Plot size in inches
        width = 10, height = 8)
      ## Add date to slide
      pptOut <- ReporteRs::addDate(doc = pptOut)
      ##Write .pptx
      ReporteRs::writeDoc(doc = pptOut, file = "pptxName")
    } else
    {message("Package ReporteRs needs to be installed to be able to export to .pptx")}
  }
}
