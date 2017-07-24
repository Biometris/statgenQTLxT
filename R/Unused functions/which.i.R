which.i <- function(i,p){
  m <- matrix(0,p,p)
  m[upper.tri(m,diag=T)] <- 1:(p*(p+1)/2)
return(as.numeric(which(m==i,arr.ind=T)))
}
#which.i(5,3)