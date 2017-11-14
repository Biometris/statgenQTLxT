#' Summary function for the class \code{GWAS}
#'
#' Gives a summary for an object of S3 class \code{GWAS}.
#'
#' @param object object of class \code{GWAS}
#' @param ... not used
#' @param environments a vector of strings or numeric indices indicating for which environment the
#' summary should be made. If \code{NULL} a summary is made for all environments.
#'
#' @export
summary.GWAS <- function(object, ..., environments = NULL) {
  ## Checks.
  if (!is.null(environments) && !is.character(environments) && !is.numeric(environments)) {
    stop("environments should be a character or numeric vector.\n")
  }
  if ((is.character(environments) && !all(environments %in% names(object$GWAResult))) ||
      (is.numeric(environments) && !all(environments %in% 1:length(object$GWAResult)))) {
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
    cat("\tData are available for", dplyr::n_distinct(GWAResult$snp), "SNPs.\n")
    if (!is.null(GWASInfo$MAF)) {
      cat("\t", dplyr::n_distinct(GWAResult$snp[is.na(GWAResult$pValue)]), "of them were not",
          "analyzed because their minor allele frequency is below", GWASInfo$MAF, "\n\n")
    }
    for (trait in traits) {
      cat("\tTrait:", trait, "\n\n")
      if (substr(GWASInfo$call[[1]], 4, 4) == "S" && !is.null(GWASInfo$GLSMethod) &&
          as.numeric(GWASInfo$GLSMethod) == 1) {
        ## Print mixed model info.
        cat("\t\tMixed model with only polygenic effects, and no marker effects:\n")
        cat("\t\tGenetic variance:", GWASInfo$varComp[[environment]][[trait]][1], "\n")
        cat("\t\tResidual variance:", GWASInfo$varComp[[environment]][[trait]][2], "\n\n")
      }
      if (!is.null(GWASInfo$thrType) && as.numeric(GWASInfo$thrType) %in% 1:3) {
        ## Print significant SNP info.
        cat("\t\tLOD-threshold:", object$thr[[environment]][trait], "\n")
        signSnpTrait <- signSnp[signSnp$trait == trait, ]
        if (!is.null(signSnpTrait)) {
          nSignSnp <- nrow(signSnpTrait[signSnpTrait$snpStatus == "significant snp", ])
          cat("\t\tNumber of significant SNPs:" , nSignSnp, "\n")
          if (nSignSnp > 0) {
            cat("\t\tSmallest p-value among the significant SNPs:",
                min(signSnpTrait[signSnpTrait$snpStatus == "significant snp", "pValue"]), "\n")
            cat("\t\tLargest p-value among the significant SNPs: ",
                max(signSnpTrait[signSnpTrait$snpStatus == "significant snp", "pValue"]),
                " (LOD-score: ", min(signSnpTrait[signSnpTrait$snpStatus == "significant snp", "LOD"]),
                ")\n\n", sep = "")
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
