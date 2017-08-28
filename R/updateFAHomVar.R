#' updateFAHomVar
#'
#' updateFAHomVar
#'
#' @inheritParams EMFA
#'
#' @param S an p x p sample covariance matrix.
#' @param m an integer. The order of the model.
#' @param maxDiag a numerical value for the maximum value of sigma2.
#'
#' @return a list containing the new value for the matrices W and P.
#'
#' @keywords internal

updateFAHomVar <- function(Y = NULL, S = NULL, m, maxDiag = 1e4) {
  if ((is.null(Y) && is.null(S)) || (!is.null(Y) && !is.null(S)))
  {stop("Either the data (Y) or the sample covariance matrix (S) must be provided.")}
  if (m != round(m) || m < 1) {stop("m needs to be a positive integer")}

  ## If S is not in imput, compute is from Y
  if (!is.null(Y)) {
    if (!is.matrix(Y)) as.matrix(Y)
    if (anyNA(Y)) {stop('Y cannot contain missing values')}
    n <- nrow(Y)
    Y <- Matrix::Matrix(scale(Y, scale = FALSE))
    S <- crossprod(Y) / n
  }
  p <- ncol(S)
  a <- eigen(S, symmetric = TRUE)

  if (m >= p) {stop("m needs to be smaller than the number of variables")}

  sigma2 <- max(mean(a$values[-(1:m)]), 1 / maxDiag)
  if (m == 1) {
    W <- matrix(a$vectors[, 1] * sqrt(a$values[1] - sigma2))
  } else {
    W <- a$vectors[, 1:m] %*% diag(sqrt(a$values[1:m] - sigma2))
  }
  return(list(W = W, P = diag(x = 1 / sigma2, nrow = p)))
}
