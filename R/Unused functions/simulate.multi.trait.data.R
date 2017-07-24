simulate.multi.trait.data <- function(K,Vg) {
  n  <- nrow(K)
  p  <- ncol(Vg)
  S  <- MatrixRoot(K)
  V  <- MatrixRoot(Vg)
  Z  <- matrix(rnorm(n=n*p),ncol=p)
  Y  <- S %*% Z %*% V
return(Y)
}
#simulate.multi.trait.data(K=K,Vg=diag(40))
