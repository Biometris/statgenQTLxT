#' Compute Square root of a symmetric, positive definite matrix
#'
#' Given an symmetric, positive definite matrix X a matrix Y is computed such that \eqn{Y^2 = X}.
#' Computation is don using eigendecomposition of X.
#'
#' @param X A symmetric, positive definite matrix
#'
#' @return A matrix Y such that \eqn{Y^2 = X}.
#'
#' @examples Z <- matrix(c(2, -1, -1, 2), nrow = 2)
#' matrixRoot(Z)

matrixRoot <- function(X) {
  if (length(dim(X)) > 2L || !is.numeric(X))
    stop("'X' must be a numeric matrix")
  if (!is.matrix(X)) X <- as.matrix(X)
  stopifnot(isSymmetric(X))
  stopifnot(is.positive.definite(X))

  XEig <- eigen(X, symmetric = TRUE)
  if (length(XEig$values) > 1) {
    XSqrt <- XEig$vectors %*% diag(sqrt(XEig$values)) %*% solve(XEig$vectors)
  } else {
    XSqrt <- XEig$vectors %*% matrix(sqrt(XEig$values)) %*% solve(XEig$vectors)
  }
  return(XSqrt)
}
