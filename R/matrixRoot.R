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
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!isSymmetric(X)) stop("X should be a symmetric matrix.")
  XEig <- eigen(X, symmetric = TRUE)
  if(any(XEig$values < 0)) stop("X should be a positive definite matrix.")
  if (length(XEig$values) > 1) {
    XSqrt <- XEig$vectors %*% diag(sqrt(XEig$values)) %*% solve(XEig$vectors)
  } else {
    XSqrt <- XEig$vectors %*% matrix(sqrt(XEig$values)) %*% solve(XEig$vectors)
  }
  return(XSqrt)
}
