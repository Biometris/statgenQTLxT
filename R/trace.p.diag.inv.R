trace.p.diag.inv <- function(x,y) {
# goal : compute directly the trace of Lambda_1^{*} and Lambda_2^{*}
# input : x,y, of type numeric
  p <- length(x)
  n <- length(y)
  z <- sapply(x,function(x_i){return(sum(1 / (rep(1,n) + x_i * y)))})
return(z)
}
#trace.p.diag.inv(1:2,1:3)
