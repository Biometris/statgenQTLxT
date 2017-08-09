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

astle <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- ncol(X)
  p <- colSums(X) / (2 * nrow(X))
  Z <- sapply(1:n, function(i) {sqrt(1 / n) * (X[, i] - 2 * p[i]) / sqrt(2 * p[i] * (1 - p[i]))})
  K <- tcrossprod(Z)
  K <- as.matrix(Matrix::nearPD(K)$mat)
}



