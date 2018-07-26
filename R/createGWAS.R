#' S3 Class GWAS
#'
#' \code{createGWAS} creates an object of S3 class GWAS containing the results
#' of a GWAS analysis.
#' \code{GWAResult} and \code{signSnp} are both optional, however at least one
#' of those should be provided as input.\cr
#' \code{summary} and \code{plot} functions are available.\cr\cr
#' \code{is.gData} tests if an \code{R} object is of class \code{gData}.
#'
#' @param GWAResult An optional data.frame or list of data.frames containing the
#' overall analysis results. Should a least contain columns \code{trait}, the
#' evaluated trait, \code{snp}, the name of the SNP, \code{chr}, the chromosome
#' number, \code{pos}, the position of the SNP on the chromosome,
#' \code{pValue}, the p-values from the analysis and \code{LOD} the LOD-score.
#' @param signSnp An optional data.frame or list of data.frames containing
#' information on the significant SNPs and optionally the SNPs close to the
#' significant SNPs. Should at least contain columns \code{trait}, the
#' evaluated trait, \code{snp}, the name of the SNP, \code{pValue}, the p-values
#' from the analysis and \code{LOD} the LOD-score.
#' @param kin An optional kinship matrix or list of chromosome specific kinship
#' matrices.
#' @param thr An optional numerical value, the threshold used in performing the
#' GWAS analysis.
#' @param GWASInfo a list containing extra information concering the GWAS
#' analysis.
#' @param x An \code{R} object
#'
#' @return \code{createGWAS} returns an object of class GWAS, a list of the
#' input items.\cr\cr
#' \code{is.GWAS} returns \code{TRUE} or \code{FALSE} depending on whether its
#' argument is a \code{GWAS} object.
#'
#' @seealso \code{\link{summary.GWAS}}, \code{\link{plot.GWAS}}
#'
#' @name GWAS
NULL

#' @rdname GWAS
#' @export
createGWAS <- function(GWAResult = NULL,
                       signSnp = NULL,
                       kin = NULL,
                       thr = NULL,
                       GWASInfo = NULL) {
  ## Check that at least one input is provided.
  if (is.null(GWAResult) && is.null(signSnp)) {
    stop("At least one of GWAResult and signSnp should be provided.\n")
  }
  ## Check GWAResults
  if (!is.null(GWAResult)) {
    if (!is.data.frame(GWAResult) &&
        !(is.list(GWAResult) &&
          all(sapply(X = GWAResult, FUN = is.data.frame)))) {
      stop("GWAResult should be a data.frame or a list data.frames.\n")
    }
    if (is.data.frame(GWAResult)) {
      ## If not a list already put data.frame in a list.
      GWAResult <- list(GWAResult)
    }
    if (!all(sapply(GWAResult, FUN = function(x) {
      all(c("trait", "snp", "chr", "pos", "pValue", "LOD") %in%
          colnames(x))}))) {
      stop(paste("GWAResult should contain columns trait, snp, chr, pos,",
                 "pValue and LOD.\n"))
    }
  }
  ## Check signSnps
  if (!all(sapply(signSnp, FUN = is.null))) {
    if (!is.data.frame(signSnp) &&
        !(is.list(signSnp) && all(sapply(X = signSnp, FUN = function(x) {
          is.null(x) || is.data.frame(x)
        })))) {
      stop("signSnp should be a data.frame or a list of data.frames.\n")
    }
    if (is.data.frame(signSnp)) {
      ## If not a list already put data.frame in a list.
      signSnp <- list(signSnp)
    }
    if (!all(sapply(signSnp, FUN = function(x) {
      is.null(x) || all(c("trait", "snp", "snpStatus", "pValue", "LOD") %in%
                        colnames(x))}))) {
      stop(paste("signSnp should contain columns trait, snp, snpStatus,",
                 "pValue and LOD.\n"))
    }
  }
  ## Check kin
  if (!is.null(kin)) {
    if (!(inherits(kin, "Matrix") || is.matrix(kin)) &&
        !(is.list(kin) && all(sapply(X = kin, FUN = function(x) {
          inherits(x, "Matrix") || is.matrix(x)})))) {
      stop("kin should be a matrix or a list of matrices.\n")
    }
  }
  ## Create GWAS object.
  GWAS <- structure(list(GWAResult = GWAResult,
                         signSnp = signSnp,
                         kinship = kin,
                         thr = thr,
                         GWASInfo = GWASInfo),
                    class = "GWAS")
  return(GWAS)
}

#' @rdname GWAS
#' @export
is.GWAS <- function(x) {
  inherits(x, "GWAS")
}

#' Summary function for the class \code{GWAS}
#'
#' Gives a summary for an object of S3 class \code{GWAS}.
#'
#' @param object An object of class \code{GWAS}
#' @param ... Not used
#' @param environments A vector of strings or numeric indices indicating for
#' which environment the summary should be made. If \code{NULL} a summary is
#' made for all environments.
#'
#' @export
summary.GWAS <- function(object, ..., environments = NULL) {
  ## Checks.
  if (!is.null(environments) && !is.character(environments) &&
      !is.numeric(environments)) {
    stop("environments should be a character or numeric vector.\n")
  }
  if ((is.character(environments) &&
       !all(environments %in% names(object$GWAResult))) ||
      (is.numeric(environments) &&
       !all(environments %in% 1:length(object$GWAResult)))) {
    stop("all environments should be in object.\n")
  }
  ## Convert character input to numeric.
  if (is.character(environments)) {
    environments <- which(names(object$GWAResult) == environments)
  }
  ## If NULL then summary of all environments.
  if (is.null(environments)) {
    environments <- 1:length(object$GWAResult)
  }
  for (environment in environments) {
    GWAResult <- object$GWAResult[[environment]]
    signSnp <- object$signSnp[[environment]]
    GWASInfo <- object$GWASInfo
    traits <- unique(GWAResult$trait)
    ## Print environment.
    cat(names(object$GWAResult)[environment], ":\n", sep = "")
    ## Print traits.
    cat("\tTraits analysed:", paste(traits, collapse = ", "), "\n\n")
    ## Print SNP numbers.
    cat("\tData are available for", length(unique(GWAResult$snp)),
        "SNPs.\n")
    if (!is.null(GWASInfo$MAF)) {
      cat("\t", length(unique(GWAResult$snp[is.na(GWAResult$pValue)])),
          "of them were not analyzed because their minor allele frequency is",
          "below", GWASInfo$MAF, "\n\n")
    }
    for (trait in traits) {
      cat("\tTrait:", trait, "\n\n")
      if (substr(GWASInfo$call[[1]], 4, 4) == "S" &&
          !is.null(GWASInfo$GLSMethod) && as.numeric(GWASInfo$GLSMethod) == 1) {
        ## Print mixed model info.
        cat("\t\tMixed model with only polygenic effects,",
            "and no marker effects:\n")
        cat("\t\tGenetic variance:",
            GWASInfo$varComp[[environment]][[trait]][1], "\n")
        cat("\t\tResidual variance:",
            GWASInfo$varComp[[environment]][[trait]][2], "\n\n")
      }
      if (!is.null(GWASInfo$thrType) && as.numeric(GWASInfo$thrType) %in% 1:3) {
        ## Print significant SNP info.
        cat("\t\tLOD-threshold:", object$thr[[environment]][trait], "\n")
        signSnpTrait <- signSnp[signSnp$trait == trait, ]
        if (!is.null(signSnpTrait)) {
          nSignSnp <-
            nrow(signSnpTrait[signSnpTrait$snpStatus == "significant snp", ])
          cat("\t\tNumber of significant SNPs:" , nSignSnp, "\n")
          if (nSignSnp > 0) {
            cat("\t\tSmallest p-value among the significant SNPs:",
                min(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "pValue"]), "\n")
            cat("\t\tLargest p-value among the significant SNPs: ",
                max(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "pValue"]),
                " (LOD-score: ",
                min(signSnpTrait[signSnpTrait$snpStatus == "significant snp",
                                 "LOD"]), ")\n\n", sep = "")
          } else {
            cat("\n")
          }
        } else {
          cat("\t\tNo significant SNPs found.","\n\n")
        }
      }
      if (!is.null(GWASInfo$genomicControl) && GWASInfo$genomicControl) {
        ## Print genomic control.
        cat("\t\tGenomic control correction was applied\n")
      } else {
        cat("\t\tNo Genomic control correction was applied\n")
      }
      if (!is.null(GWASInfo$inflationFactor)) {
        cat("\t\tGenomic control inflation-factor:",
            round(GWASInfo$inflationFactor[[environment]][trait], 3), "\n\n")
      }
    }
  }
}

#' Plot function for the class \code{GWAS}
#'
#' Creates a plot of an object of S3 class \code{GWAS}. Three types of plots
#' can be made:
#' \itemize{
#' \item{a manhattan plot using the function \code{\link{manhattanPlot}}}
#' \item{a qq plot using the function \code{\link{qqPlot}}}
#' \item{a qtl plot using the function \code{\link{qtlPlot}}}
#' }
#'
#' When making a manhattan plot all markers in the \code{GWAResult} data.frame
#' in the \code{GWAS} object are plotted. Significant SNPs are taken from the
#' \code{signSnp} data.frame in the \code{GWAS} object and marked in the plot.
#' Also the LOD-threshold is taken from the \code{GWAS} object and plotted as
#' a horizontal line. Both \code{signSnp} and \code{thr} can be left empty and
#' will then be ignored in the plot.\cr
#' \code{...} can be used to pass extra arguments to the actual plotting
#' functions. See those fuctions for details.
#'
#' @param x An object of class \code{GWAS}.
#' @param ... further arguments to be passed on to the actual plotting
#' functions.
#' @param type A character string indicating the type of plot to be made. One
#' of "manhattan", "qq" and "qtl".
#' @param environment A character string or numeric index indicating for which
#' environment the plot should be made. If \code{x} only contains results for
#' one trait \code{environment} may be \code{NULL}.
#' @param trait A character string indicating for which trait the results
#' should be plotted. For \code{type} "qtl" all traits are plotted. If \code{x}
#' only contains results for one trait \code{trait} may be \code{NULL}.
#'
#' @seealso \code{\link{manhattanPlot}}, \code{\link{qqPlot}},
#' \code{\link{qtlPlot}}
#'
#' @export
plot.GWAS <- function(x, ...,
                      type = c("manhattan", "qq", "qtl"),
                      environment = NULL,
                      trait = NULL) {
  type <- match.arg(type)
  dotArgs <- list(...)
  ## Checks.
  if (!is.null(environment) && !is.character(environment) &&
      !is.numeric(environment)) {
    stop("Environment should be a character or numerical value.\n")
  }
  if ((is.character(environment) && !environment %in% names(x$GWAResult)) ||
      (is.numeric(environment) && !environment %in% 1:length(x$GWAResult))) {
    stop("Environment should be in x.\n")
  }
  ## Convert character input to numeric.
  if (is.character(environment)) {
    environment <- which(names(x$GWAResult) == environment)
  }
  ## If NULL then summary of all environment.
  if (is.null(environment)) {
    if (length(x$GWAResult) != 1) {
      stop(paste("Environment not supplied but multiple environments",
                 "detected in data.\n"))
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
          stop("Trait not supplied but multiple traits detected in data.\n")
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
    ## Select specific chromosome(s) for plotting.
    if (!is.null(dotArgs$chr)) {
      GWAResult <- GWAResult[GWAResult$chr %in% dotArgs$chr, ]
      if (nrow(GWAResult) == 0) {
        stop("Select at least one valid chromosome for plotting.\n")
      }
    }
    ## Select markers with sufficiently high lod for plotting.
    if (!is.null(dotArgs$lod)) {
      GWAResult <- GWAResult[GWAResult$LOD > dotArgs$lod, ]
      if (nrow(GWAResult) == 0) {
        stop(paste("No chromosomes selected for plotting. Please check",
                   "value of lod.\n"))
      }
    }
    chrBnd <- aggregate(x = GWAResult$pos, by = list(GWAResult$chr), FUN = max)
    ## Compute cumulative positions.
    addPos <- data.frame(chr = chrBnd[, 1],
                         add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                         stringsAsFactors = FALSE)
    map <- GWAResult[, c("snp", "chr", "pos", "LOD")]
    map <- merge(map, addPos, by = "chr")
    map$cumPos <- map$pos + map$add
    ## Extract numbers of significant SNPs.
    if (!is.null(signSnp)) {
      if (is.null(dotArgs$yThr)) {
        signSnpNr <- which(map$snp %in%
                             signSnp$snp[as.numeric(signSnp$snpStatus) == 1])
      } else {
        signSnpNr <- which(map$LOD > dotArgs$yThr)
      }
    } else {
      signSnpNr <- integer()
    }
    if (!is.null(dotArgs$plotType)) {
      plotType = dotArgs$plotType
    } else {
      plotType = "p"
    }
    if (is.null(dotArgs$yThr)) {
      if (is.null(x$thr[[environment]][trait])) {
        yThr <- Inf
      } else {
        yThr <- x$thr[[environment]][trait]
      }
    } else {
      yThr <- dotArgs$yThr
    }
    ## Create manhattan plot.
    do.call(manhattanPlot,
            args = c(list(xValues = map$cumPos, yValues = map$LOD,
                          map = map[, -which(colnames(map) == "LOD")],
                          plotType = plotType, xSig = signSnpNr,
                          chrBoundaries = chrBnd[, 2], yThr = yThr),
                     dotArgs[!(names(dotArgs) %in% c("plotType", "yThr",
                                                     "lod", "chr"))]
            ))
  } else if (type == "qq") {
    ## Create qq-plot
    qqPlot(pValues = na.omit(GWAResult$pValue), ...)
  } else if (type == "qtl") {
    if (is.null(signSnp)) {
      stop("No significant SNPs in signSnp. No plot can be made.\n")
    }
    qtlPlot(data = signSnp,
            map = GWAResult[!is.na(GWAResult$pos), c("chr", "pos")],
            ...)
  }
}
