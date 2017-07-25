# X = n x p matrix of marker scores
GRM <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- ncol(X)
  K <- tcrossprod(scale(X)) / n
  K <- as.matrix(Matrix::nearPD(K)$mat)
  return(K)
}

