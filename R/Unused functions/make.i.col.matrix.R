make.i.col.matrix <- function(i,p) {
  i.matrix        <- matrix(0,p,1)
  i.matrix[i,1]  <- 1
return(i.matrix)
}
# make.i.col.matrix(3,p=7)
