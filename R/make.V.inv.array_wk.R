make.V.inv.array <- function(Vg, Ve, Dk) {
  n <- length(Dk)
  p <- ncol(Vg)
  V.inv.array <- array(dim = c(n, p, p))
  for (i in 1:n) {
    V.inv.array[i, , ] <- solve(Dk[i] * Vg + Ve)
  }
  return(V.inv.array)
}
