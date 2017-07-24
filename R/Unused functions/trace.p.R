trace.p <- function(np.matrix,n,p) {
  result <- matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      result[i,j] <- sum(diag(np.matrix[(i-1)*p + 1:p,(j-1)*p + 1:p]))
    }
  }
return(result)
}