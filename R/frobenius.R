frobenius <- function(M) {
  return(sum(diag(t(M) %*% M)))
}