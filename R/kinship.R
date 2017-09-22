#' Functions for calculating kinship matrices
#'
#' A collection of functions for calculating kinship matrices using different algorithms. The following
#' algorithms are included: astle (Astle and Balding, 2009), GRM, Identity By State (IBS) and
#' VanRaden (VanRaden, 2008).
#'
#' In all algorithms the input matrix \code{X} is first cleaned, i.e. markers with a variance of 0 are
#' excluded from the calculation of the kinship matrix. Then some form of scaling is done which differs
#' per algorithm. This gives a scaled matrix \code{Z}. The matrix \eqn{ZZ^t / denominator} is returned.
#' By default the denominator is equal to the number of columns in \code{Z} for \code{astle}, \code{GRM}
#' and \code{IBS} and \eqn{2* p * (1-p)} where \eqn{p = colSums(X) / (2 * nrow(X))} for \code{vanRaden}.
#' This denominator can be overwritten by the user, e.g. when computing kinship matrices by splitting
#' \code{X} in smaller matrices and then adding the results together in the end.
#'
#' @param X an n x m marker matrix with genotypes in the rows and markers in the columns.
#' @param denominator an numeric value. See details.
#'
#' @return an n x n kinship matrix.
#'
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic relatedness in genetic
#' association studies, Stat. Sci., November 2009, Vol. 24, no. 4, p. 451–471.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic predictions. J Dairy Sci,
#' November 2008, Vol. 91 p. 4414–4423.
#'
#' @name kinship
NULL

#' @examples X <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1), nrow = 3)
#' astle(X)
#' GRM(X)
#' IBS(X)
#' vanRaden(X)

#' @rdname kinship
#' @export
astle <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  ## Remove markers with variance 0.
  X <- X[, apply(X, 2, var) != 0, drop = FALSE]
  ## Scale X.
  p <- colSums(X) / (2 * nrow(X))
  Z <- scale(X, center = 2 * p, scale = sqrt(2 * p * (1 - p)))
  ## Compute denominator.
  if (is.null(denominator)) denominator <- ncol(Z)
  return(Matrix::tcrossprod(as(Z, "dgeMatrix")) / denominator)
}

#' @rdname kinship
#' @export
GRM <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  ## Remove markers with variance 0.
  X <- X[, apply(X, 2, var) != 0, drop = FALSE]
  ## Scale X.
  Z <- scale(X)
  ## Compute denominator.
  if (is.null(denominator)) denominator <- ncol(Z)
  return(Matrix::tcrossprod(as(Z, "dgeMatrix")) / denominator)
}

#' @rdname kinship
#' @export
IBS <- function(X,
  denominator = NULL) {
  ## Remove markers with variance 0.
  X <- X[, apply(X, 2, var) != 0, drop = FALSE]
  ## Compute denominator.
  if (is.null(denominator)) denominator <- ncol(X)
  return((Matrix::tcrossprod(X) + Matrix::tcrossprod(1 - X)) / denominator)
}

#' @rdname kinship
#' @export
vanRaden <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  ## Remove markers with variance 0.
  X <- X[, apply(X, 2, var) != 0, drop = FALSE]
  ## Scale X.
  p <- colSums(X) / (2 * nrow(X))
  Z <- scale(X, center = 2 * p, scale = FALSE)
  ## Compute denominator.
  if (is.null(denominator)) denominator <- 2 * sum(p * (1 - p))
  return(Matrix::tcrossprod(as(Z, "dgeMatrix")) / denominator)
}


