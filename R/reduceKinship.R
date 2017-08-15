#' Reduce the kinship matrix
#'
#' The kinship matrix is reduced using nPca eigenvectors of K.
#'
#' @inheritParams runMultiTraitGwas
#'
#' @param nPca an integer, the number of eigenvectors used for reducing the kinship matrix.
#'
#' @return the reduced kinship matrix

reduceKinship <- function(K, nPca) {
  w <- eigen(K, symmetric = TRUE)
  U <- w$vectors[, 1:nPca]
  S <- diag(w$values[1:nPca])
  KRed <- U %*% S %*% t(U)
  return(KRed)
}
