reduce_kinship <- function(K,n) {
      ev <- eigen(K)
      U  <- ev$vectors[,1:n]
      S  <- diag(ev$values[1:n])
      K.red <- U %*% S %*% t(U)
  return(K.red)
}