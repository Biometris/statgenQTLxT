#' Summary function for the class \code{GData}
#'
#' Gives a summary for an object of S3 class \code{GData}.
#'
#' @param object An object of class \code{GData}
#' @param ... Not used
#'
#' @export
summary.gData <- function(object, ...) {
  map <- object$map
  markers <- object$markers
  pheno <-  object$pheno
  kinship <- object$kinship
  covar <- object$covar
  if (!is.null(map)) {
    cat("map\n")
    cat("\tNumber of markers:", nrow(map), "\n")
    cat("\tNumber of chromosomes:", dplyr::n_distinct(map$chr), "\n\n")
  }
  if (!is.null(markers)) {
    cat("markers\n")
    cat("\tNumber of markers:", ncol(markers), "\n")
    cat("\tNumber of genotypes:", nrow(markers), "\n")
    cat("\tContent:\n")
    tab <- (round(prop.table(table(as.vector(markers), useNA = "always")), 2))
    cat("\t", names(tab), "\n")
    cat("\t", tab, "\n\n")
  }
  if (!is.null(pheno)) {
    cat("pheno\n")
    cat("\tNumber of environments:", length(pheno), "\n\n")
    for (i in 1:length(pheno)) {
      if (!is.null(names(pheno)[i])) {
        cat("\t", names(pheno)[i], ":\n", sep = "")
      } else {
        cat("\tEnvironment ", i, ":\n", sep = "")
      }
      cat("\t\tNumber of traits:", ncol(pheno[[i]]) - 1, "\n")
      cat("\t\tNumber of genotypes:", dplyr::n_distinct(pheno[[i]]$genotype),
          "\n\n")
      print(summary(pheno[[i]][, -1]))
      cat("\n")
    }
  }
  if (!is.null(covar)) {
    cat("covar\n")
    cat("\tNumber of covariates:", ncol(covar), "\n")
    print(summary(covar))
  }
}
