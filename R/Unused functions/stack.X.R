stack.X <- function(X,p=1) {
  require(Matrix)
  r <- list()
  for (i in 1:p) {
    r[[i]] <- X
  }
return(as.matrix(bdiag(r)))
}