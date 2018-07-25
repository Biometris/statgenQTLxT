#' Update W and P in EMFA algorithm for homogeneous variance.
#'
#' Updata W and P used in the iteration process in the EMFA algorithm in case
#' the variance is homogeneous.
#'
#' @inheritParams EMFA
#'
#' @param S An p x p sample covariance matrix.
#' @param m An integer. The order of the model.
#' @param maxDiag A numerical value for the maximum value of sigma2.
#'
#' @return A list containing the new value for the matrices W and P.
#'
#' @keywords internal

updateFAHomVar <- function(Y = NULL,
                           S = NULL,
                           m,
                           maxDiag = 1e4) {
  if ((is.null(Y) && is.null(S)) || (!is.null(Y) && !is.null(S))) {
    stop(paste("Either the data (Y) or the sample covariance matrix (S)",
               "must be provided.\n"))
  }
  if (m != round(m) || m < 1) {
    stop("m needs to be a positive integer")
  }
  ## If S is not in imput, compute S from Y.
  if (!is.null(Y)) {
    if (anyNA(Y)) {
      stop("Y cannot contain missing values.\n")
    }
    n <- nrow(Y)
    Y <- as(scale(Y, scale = FALSE), "dgeMatrix")
    S <- Matrix::crossprod(Y) / n
  }
  p <- ncol(S)
  a <- eigen(S, symmetric = TRUE)
  if (m >= p) {
    stop("m needs to be smaller than the number of variables.\n")
  }
  sigma2 <- max(mean(a$values[-(1:m)]), 1 / maxDiag)
  ## Split cases for extra robustness.
  if (m == 1) {
    W <- as(a$vectors[, 1] * sqrt(a$values[1] - sigma2), "dgeMatrix")
  } else {
    W <- a$vectors[, 1:m] %*% Matrix::Diagonal(x = sqrt(a$values[1:m] - sigma2))
  }
  return(list(W = W, P = Matrix::Diagonal(n = p, x = 1 / sigma2)))
}
