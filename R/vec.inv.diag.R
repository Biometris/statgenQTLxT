vec.inv.diag <- function(x,y) {
# goal : compute directly the trace of Lambda_1^{*} and Lambda_2^{*}
# input : x,y, of type numeric
# x=1:2;y=1:3
  p <- length(x)
  n <- length(y)
  z <- sapply(x,function(x_i){return(1 / (rep(1,n) + x_i * y))})
return(z)
}
# vec.inv.diag(1:2,1:3)
