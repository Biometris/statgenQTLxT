# Taken from : http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
MatrixRoot <- function(x) { # assumes that x is symmetric
  x.eig <- eigen(as.matrix(x),symmetric=TRUE)
  if (length(x.eig$values) > 1) {
    x.sqrt <- x.eig$vectors %*% diag(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  } else {
    x.sqrt <- x.eig$vectors %*% matrix(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  }
  return(x.sqrt)
}