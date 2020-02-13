#' Fast Generalized Least Squares algoritm for IBD based QTL mapping
#'
#' Compute statistics for the Generalized Least Squares (GLS) F-test for
#' IBD based QTL Mapping.
#'
#' @param y A numeric vector of length n of phenotypic scores. No missing
#' values allowed.
#' @param X An array of marker probabilities.
#' @param Sigma An n x n covariance matrix. No missing values allowed.
#' @param covs An n x c matrix of covariates (NOT including an intercept).
#' No missing values allowed.
#' @param ref An integer indicating the allele to use as reference allele in
#' the computations.
#'
#' @keywords internal
fastGLSIBD <- function(y,
                       X,
                       Sigma,
                       covs = NULL,
                       ref,
                       nCores = NULL) {
  ## Check class and missing values.
  if (missing(y) || !(inherits(y, "matrix") || is.numeric(y)) || anyNA(y)) {
    stop("y should be a numeric vector without missing values.\n")
  }
  if (missing(Sigma) || !(inherits(Sigma, "matrix") || is.matrix(Sigma)) ||
      anyNA(Sigma)) {
    stop("Sigma should be a matrix without missing values.\n")
  }
  if (!is.null(covs) && (!(inherits(covs, "matrix") || is.matrix(covs)) ||
                         anyNA(covs))) {
    stop("covs should be a numeric vector without missing values.\n")
  }
  n <- length(y)
  ## Check dimensions.
  if (nrow(Sigma) != n || ncol(Sigma) != n) {
    stop("The number of elements in y should be identical to the",
         "number of rows and columns in Sigma.\n")
  }
  if (!is.null(covs) && nrow(covs) != n) {
    stop("The number of elements in y should be identical to the",
         "number of rows in covs.\n")
  }
  resCpp <- fastGLSIBDCPP(X, y, Sigma, ref, covs, nCores = nCores)
  GLS <- cbind(resCpp$pVal, resCpp$RLR2, t(resCpp$beta2))
  rownames(GLS) <- colnames(X)
  colnames(GLS) <- c("pValue", "RLR2", dimnames(X)[[3]][-ref])
  return(GLS)
}
