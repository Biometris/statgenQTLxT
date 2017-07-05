make.ij.matrix <- function(i,j,p,sym=TRUE) {
  ij.matrix        <- matrix(0,p,p)
  ij.matrix[i,j]   <- ij.matrix[i,j] + 1
  if (sym) {
    ij.matrix[j,i]   <- ij.matrix[j,i] + 1
  }
return(ij.matrix)
}