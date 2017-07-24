matrix.to.par <- function(M) {
  stopifnot(ncol(M)==nrow(M))
return(matrix(M[upper.tri(M,diag=T)]))
}