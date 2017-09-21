#' Compute Square root of a symmetric, positive definite matrix
#'
#' Given an symmetric, positive definite matrix X a matrix Y is computed such that \eqn{Y^2 = X}.
#' Computation is done using eigendecomposition of X.
#'
#' @param X a symmetric, positive definite matrix.
#'
#' @return a matrix Y such that \eqn{Y^2 = X}.
#'
#' @keywords internal

matrixRoot <- function(X) {
  #if (!is.matrix(X)) X <- as.matrix(X)
  if (!Matrix::isSymmetric(X)) stop("X should be a symmetric matrix.")
  XEig <- eigen(X, symmetric = TRUE)
  if(any(XEig$values < 0)) stop("X should be a positive definite matrix.")
  if (length(XEig$values) > 1) {
    XSqrt <- as(XEig$vectors, "dgeMatrix") %*% Matrix::Diagonal(x = sqrt(XEig$values)) %*%
      Matrix::solve(XEig$vectors)
  } else {
    XSqrt <-  as(XEig$vectors, "dgeMatrix") %*% Matrix::Matrix(sqrt(XEig$values)) %*%
      Matrix::solve(XEig$vectors)
  }
  return(XSqrt)
}
