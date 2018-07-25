#' Compute Square root of a symmetric, positive definite matrix
#'
#' Given a symmetric, positive definite matrix X a matrix Y is computed such
#' that \eqn{Y^2 = X}. Computation is done using eigendecomposition of X.
#'
#' @param X A symmetric, positive definite matrix.
#'
#' @return A matrix Y such that \eqn{Y^2 = X}.
#'
#' @keywords internal

matrixRoot <- function(X) {
  #if (!is.matrix(X)) X <- as.matrix(X)
  if (!Matrix::isSymmetric(X)) {
    stop("X should be a symmetric matrix.\n")
  }
  XEig <- eigen(X, symmetric = TRUE)
  if (any(XEig$values < 0)) {
    stop("X should be a positive definite matrix.\n")
  }
  if (length(XEig$values) > 1) {
    XSqrt <- as(XEig$vectors, "dgeMatrix") %*%
      Matrix::Diagonal(x = sqrt(XEig$values)) %*% Matrix::solve(XEig$vectors)
  } else {
    XSqrt <- as(XEig$vectors, "dgeMatrix") %*%
      Matrix::Matrix(sqrt(XEig$values)) %*% Matrix::solve(XEig$vectors)
  }
  return(XSqrt)
}
