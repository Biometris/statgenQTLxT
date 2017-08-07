# X = n x p matrix of marker scores
GRM <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  X <- scale(X)
  X <- X[, !is.na(colSums(X))]
  n <- ncol(X)
  K <- tcrossprod(X) / n
  K <- as.matrix(Matrix::nearPD(K)$mat)
  return(K)
}


