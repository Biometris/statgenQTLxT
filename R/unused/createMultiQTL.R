#' Plot function for the class \code{GWAS}
#'
#' Creates a plot of an object of S3 class \code{GWAS}. The following types of
#' plot can be made:
#' \itemize{
#' \item{a manhattan plot, i.e. a plot of LOD-scores per SNP}
#' \item{a qq plot of observed LOD-scores versus expected LOD-scores}
#' \item{a qtl plot of effect sizes and directions for multiple traits (not
#' for IBD based QTL mapping)}
#' \item{a matrix plot of effect sizes and directions per founder (only for IBD
#' based QTL mapping)}
#' }
#' Manhattan plots, qq plots and matrix plots are made for a single trait which
#' should be indicated using the parameter \code{trait} unless the analysis was
#' done for only one trait in which case it is detected automatically. The qtl
#' plot will plot all traits analysed.\cr
#' See details for a detailed description of the plots and the plot options
#' specific to the different plots.
#'
#' @section Manhattan Plot:
#' A LOD-profile of all marker positions and corresponding LOD-scores is
#' plotted. Significant markers are highlighted with red dots. By default these
#' are taken from the result of the GWAS analysis however the LOD-threshold for
#' significant parameters may be modified using the parameter \code{yThr}. The
#' treshold is plotted as a horizontal line. If there are previously known
#' marker effect, false positives and true negatives can also be marked.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{-log10(p)}}
#' \item{\code{effects}}{A character vector, indicating which SNPs correspond
#' to a real (known) effect. Used for determining true/false positives and
#' false negatives. True positives are colored green, false positives orange and
#' false negatives yellow.}
#' \item{\code{colPalette}}{A color palette used for plotting. Default
#' coloring is done by chromosome, using black and grey.}
#' \item{\code{yThr}}{A numerical value for the LOD-threshold. The value from
#' the GWAS analysis is used as default.}
#' \item{\code{signLwd}}{A numerical value giving the thickness of the
#' points that are false/true positives/negatives. Default = 0.6}
#' \item{\code{lod}}{A positive numerical value. For the SNPs with a LOD-value
#' below this value, only 5\% is plotted. The chance of a SNP being plotting is
#' proportional to its LOD-score. This option can be useful when plotting a
#' large number of SNPs.}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default, all
#' chromosomes are plotted. Using this option allows restricting the plot to a
#' subset of chromosomes.}
#' }
#'
#' @section QQ Plot:
#' From the LOD-scores calculated in the GWAS analysis, a qq-plot is generated with
#' observed LOD-scores versus expected LOD-scores. Code is adapted from
#' Segura et al. (2012).
#'
#' @section QTL Plot:
#' A plot of effect sizes for the significant SNPs found in the GWAS analysis
#' is created. Each horizontal line contains QTLs of one trait,
#' phenotypic trait or trial. Optionally, vertical white lines can indicate
#' chromosome subdivision, genes of interest, known QTL, etc. Circle diameters
#' are proportional to the absolute value of allelic effect. Colors indicate the
#' direction of the effect: green when the allele increases the trait value,
#' and blue when it decreases the value.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{normalize}}{Should the snpEffect be normalized? Default =
#' \code{FALSE}}
#' \item{\code{sortData}}{Should the data be sorted before plotting? Either
#' \code{FALSE}, if no sorting should be done, or a character string indicating
#' the data column to use for sorting. Default =
#' \code{FALSE}}
#' \item{\code{binPositions}}{An optional data.frame containing at leasts two
#' columns, chr(omosome) and pos(ition). Vertical lines are plotted at those
#' positions. Default = \code{NULL}}
#' \item{\code{printVertGrid}}{Should default vertical grid lines be plotted.
#' Default = \code{TRUE}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{"Traits"}}
#' \item{\code{yThr}}{A numerical value for the LOD-threshold. The value from
#' the GWAS analysis is used as default.}
#' \item{\code{chr}}{A vector of chromosomes to be plotted. By default all
#' chromosomes are plotted. Using this option this can be restricted to a
#' subset of chromosomes.}
#' \item{\code{exportPptx}}{Should the plot be exported to a .pptx file?
#' Default = \code{FALSE}}
#' \item{\code{pptxName}}{A character string, the name of the .pptx file to
#' which the plot is exported. Ignored if exportPptx = \code{FALSE}.}
#' }
#'
#' @section Matrix Plot:
#' A plot of effect sizes for each of the founders for the significant SNPs
#' found in the IBD based QTL mapping is created.\cr
#' Extra parameter options:
#' \describe{
#' \item{\code{xLab}}{A character string, the x-axis label. Default =
#' \code{"Chromosomes"}}
#' \item{\code{yLab}}{A character string, the y-axis label. Default =
#' \code{"Parents"}}
#' }
#'
#' @param x An object of class \code{SIM}.
#' @param ... further arguments to be passed on to the actual plotting
#' functions.
#' @param plotType A character string indicating the type of plot to be made.
#' One of "manhattan", "qq","qtl" and "matrix". "qtl" plots cannot be made for
#' results of IBD based QTL mapping. "matrix" plot can only be made for IBD based QTL
#' mapping.
#' @param trial A character string or numeric index indicating for which
#' trial the plot should be made. If \code{x} only contains results for
#' one trial, \code{trial} may be \code{NULL}.
#' @param trait A character string indicating for which trait the results
#' should be plotted. For \code{type} "qtl" all traits are plotted. If \code{x}
#' only contains results for one trait, \code{trait} may be \code{NULL}.
#' @param reference A named numerical vector used as reference value for
#' plotting effect size. The default value for each trait is the mean raw trait
#' value (available as attribute in the multiQTL object).
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a list of ggplot objects is invisibly returned.
#'
#' @references Millet et al. (2016) Genome-wide analysis of yield in Europe:
#' Allelic effects vary with drought and heat scenarios. Plant Physiology,
#' October 2016, Vol. 172, p. 749–764
#' @references Segura V, Vilhjálmsson BJ, Platt A, et al. An efficient
#' multi-locus mixed model approach for genome-wide association studies in
#' structured populations. Nature genetics. 2012;44(7):825-830.
#' doi:10.1038/ng.2314.
#'
#' @export
plot.multiQTL <- function(x,
                          ...,
                          plotType = c("manhattan", "qtlEff", "trtEff", "matrix"),
                          trial = NULL,
                          trait = NULL,
                          reference = NULL,
                          output = TRUE) {
  plotType <- match.arg(plotType)
  if (plotType %in% c("manhattan")) {
    statgenGWAS:::plot.GWAS(x = x, ... = ..., plotType = plotType,
                            trial = trial, trait = trait, type = "lines",
                            output = output)
  } else {
    dotArgs <- list(...)
    ## Get results.
    GWAResult <- x$GWAResult[[1]]
    ## Get peaks.
    QTLS <- x$peaks
    ## Get mean effect per trait.
    trtMeans <- attr(x, "trtMeans")
    ## Restrict reference to traits in trtMeans.
    reference <- reference[names(reference)  %in% names(trtMeans)]
    if (length(reference) > 0) {
      ## Replace by reference if specified.
      trtMeans[names(reference)] <- reference
    }
    QTLS[["trtMean"]] <- trtMeans[match(QTLS[["trait"]], names(trtMeans))]
    QTLS[["effect"]] <- QTLS[["effect"]] / QTLS[["trtMean"]]
    QTLS[["effectSe"]] <- QTLS[["effectSe"]] / QTLS[["trtMean"]]
    ## Convert snp to factor to assure order matches the order on genome.
    QTLS[["snp"]] <- factor(QTLS[["snp"]], levels = unique(QTLS[["snp"]]))
    if (plotType == "qtlEff") {
      ## Create plot.
      ggplot2::ggplot(data = QTLS,
                      mapping = ggplot2::aes_string(x = "snp", y = "effect",
                                                    fill = "trait")) +
        ggplot2::geom_bar(stat = "identity",
                          position = ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "effect - effectSe",
                                                   ymax = "effect + effectSe"),
                               position = ggplot2::position_dodge()) +
        ggplot2::labs(y = "effect as % of reference")

    } else if (plotType == "trtEff") {
      ## Create plot.
      ggplot2::ggplot(data = QTLS,
                      mapping = ggplot2::aes_string(x = "trait", y = "effect",
                                                    fill = "snp")) +
        ggplot2::geom_bar(stat = "identity",
                          position = ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "effect - effectSe",
                                                   ymax = "effect + effectSe"),
                               position = ggplot2::position_dodge()) +
        ggplot2::labs(y = "effect as % of reference")
    } else if (plotType == "matrix") {
      signSnp <- GWAResult
      signSnp[["qtl"]] <- signSnp[["snp"]] %in% QTLS[["snp"]]
      signSnp[["sign"]] <- abs(signSnp[["effect"]]) /
        signSnp[["effectSe"]] > 1.96
      signQTLS <- signSnp[signSnp[["qtl"]] & signSnp[["sign"]], ]
      signQTL <- lapply(X = 1:nrow(signQTLS), FUN = function(j) {
        z <- signQTLS[j, ]
        chrQTL <- signSnp[signSnp[["chr"]] == z[["chr"]] &
                            signSnp[["trait"]] == z[["trait"]], ]
        signChrQTL <- chrQTL[chrQTL[["sign"]], ]
        res <- lapply(X = 1:nrow(signChrQTL), FUN = function(i) {
          y <- signChrQTL[i, ]
          if ((y[["pos"]] <= z[["pos"]] &&
               all(chrQTL[chrQTL[["pos"]] >= y[["pos"]] &
                          chrQTL[["pos"]] <= z[["pos"]], ][["sign"]])) ||
              (y[["pos"]] > z[["pos"]] &&
               all(chrQTL[chrQTL[["pos"]] <= y[["pos"]] &
                          chrQTL[["pos"]] >= z[["pos"]], ][["sign"]]))) {
            y
          }
        })
        do.call(rbind, res)
      })
      signSnp <- do.call(rbind, signQTL)
      ## Compute chromosome boundaries.
      GWAResult <- GWAResult[!is.na(GWAResult$pos), ]
      ## Select specific chromosome(s) for plotting.
      if (!is.null(dotArgs$chr)) {
        GWAResult <- GWAResult[GWAResult$chr %in% dotArgs$chr, ]
        if (nrow(GWAResult) == 0) {
          stop("Select at least one valid chromosome for plotting.\n")
        }
      }
      GWAResComp <- GWAResult[GWAResult[["trait"]] == GWAResult[["trait"]][1], ]
      chrBnd <- aggregate(x = GWAResComp$pos, by = list(GWAResComp$chr), FUN = max)
      ## Compute cumulative positions.
      addPos <- data.frame(chr = chrBnd[, 1],
                           add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                           stringsAsFactors = FALSE)
      map <- GWAResComp[, c("snp", "chr", "pos", "LOD")]
      map <- merge(map, addPos, by = "chr")
      map$cumPos <- map$pos + map$add
      do.call(matrixPlot,
              args = c(list(effectDat = GWAResult, signSnp = signSnp,
                            map = map, chrBoundaries = chrBnd,
                            founders = x$GWASInfo$founders, output = output),
                       dotArgs[!(names(dotArgs) %in% c("effectDat", "signSnp",
                                                       "map", "chrBoundaries",
                                                       "founders"))]))
    }
  }
}
