# X = n x p matrix of marker scores
GRM <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  X <- X[, !colSums(X) %in% c(0, 1)]
  Z <- scale(X)
  if (is.null(denominator)) denominator <- ncol(Z)
  return(tcrossprod(X) / denominator)
}

astle <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  X <- X[, !colSums(X) %in% c(0, 1)]
  p <- colSums(X) / (2 * nrow(X))
  Z <- scale(X, center = 2 * p, scale = sqrt(2 * p * (1 - p)))
  if (is.null(denominator)) denominator <- ncol(Z)
  return(tcrossprod(Z) / denominator)
}

vanRaden <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  X <- X[, !colSums(X) %in% c(0, 1)]
  p <- colSums(X) / (2 * nrow(X))
  alpha <- 2 * sum(p * (1 - p))
  Z <- sqrt(1 / alpha) * scale(X, center = 2 * p, scale = FALSE)
  return(tcrossprod(Z))
}

IBS <- function(X,
  denominator = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  X <- X[, !colSums(X) %in% c(0, 1)]
  if (is.null(denominator)) denominator <- ncol(X)
  return((tcrossprod(X) + tcrossprod(1 - X)) / denominator)
}

