reduceKinship <- function(K, n) {
  w <- eigen(K, symmetric = TRUE)
  U <- w$vectors[, 1:n]
  S <- diag(w$values[1:n])
  KRed <- U %*% S %*% t(U)
  return(KRed)
}
