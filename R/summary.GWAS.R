#' Summary function for the class \code{GWAS}
#'
#' Gives a summary for an object of S3 class \code{GWAS}.
#'
#' @param object object of class \code{GWAS}
#' @param ... not used
#'
#' @export
summary.GWAS <- function(object, ...) {
  GWAResult <- object$GWAResult[[1]]
  signSnp <- object$signSnp[[1]]
  GWASInfo <- object$GWASInfo
  ## Print traits
  cat("Traits analysed:", paste(unique(GWAResult$trait), collapse = ", "), "\n\n")
  ## Print SNP numbers.
  cat("Data are available for", length(unique(GWAResult$snp)), "SNPs.\n")
  if (!is.null(GWASInfo$MAF)) {
    cat(length(unique(GWAResult$snp[is.na(GWAResult$pValue)])), "of them were not analyzed because their",
      "minor allele frequency is below", GWASInfo$MAF, "\n\n")
  }
  if (as.numeric(GWASInfo$GLSMethod) == 1) {
    ## Print mixed model info.
    cat("Mixed model with only polygenic effects, and no marker effects:\n")
    cat("Genetic variance:", GWASInfo$varComp[[1]][1], "\n")
    cat("Residual variance:", GWASInfo$varComp[[1]][2], "\n\n")
  }
  if (as.numeric(GWASInfo$thrType) %in% 1:3) {
    ## Print significant SNP info.
    cat("LOD-threshold:", object$thr, "\n")
    if (!is.null(signSnp)) {
      cat("Number of significant SNPs:" , nrow(signSnp[signSnp$snpStatus == "significant snp", ]), "\n")
      cat("Smallest p-value among the significant SNPs:",
        min(signSnp[signSnp$snpStatus == "significant snp", "pValue"]), "\n")
      cat("Largest p-value among the significant SNPs: ",
        max(signSnp[signSnp$snpStatus == "significant snp", "pValue"]),
        " (LOD-score: ", min(signSnp[signSnp$snpStatus == "significant snp", "LOD"]), ")\n", sep = "")
    } else {
      cat("No significant SNPs found.","\n")
    }
    if (GWASInfo$genomicControl) {
      ## Print genomic control.
      cat("\nGenomic control correction was applied\n")
      cat("Genomic control inflation-factor:",
        paste(round(GWASInfo$inflationFactor, 2), collapse = ", "), "\n\n")
    } else {
      cat("\nNo Genomic control correction was applied\n")
    }
  }
}
