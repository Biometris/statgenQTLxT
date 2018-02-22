make.V.inv.array <- function(Vg,Ve,Dk) {

  stopifnot(ncol(Vg)==nrow(Vg))
  stopifnot(ncol(Ve)==nrow(Ve))
  stopifnot(ncol(Vg)==ncol(Ve))

  n <- ncol(Dk)
  p <- ncol(Vg)

  V.inv.array <- array(0,dim=c(n,p,p))
  for (i in 1:n) {
    V.inv.array[i,,] <- solve(diag(Dk)[i] * Vg + Ve)
  }
return(V.inv.array)
}