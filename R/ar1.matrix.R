ar1.matrix <- function(n=5,rho=0.5) {
  H <- abs(outer(1:n,1:n, "-"))
  V <- rho^H
return(V)
}