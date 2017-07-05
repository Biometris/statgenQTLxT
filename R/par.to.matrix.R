# Assuming p traits, the genetic- and environmental covariance matrices
# are both p x p, each having  Q = p*(p+1)/2 parameters. For a symmetric
# p x p matrix A with elements a_{ij}, we define the corresponding parameter
# vector to be (a_{11},a_{12},a_{22},a_{13},a_{23},a_{33},...,a_{pp})
# i.e. the upper triangular part, including the diagonal; column by column
# The function matrix.to.par (given below) converts a p x p matrix to the
# Q x 1 column vector containing the parameters
# The function par.to.matrix performs the conversion in the other direction

par.to.matrix <- function(v) {
  # v is a vector of type numeric, or a column/row vector of type matrix
  v <- as.numeric(v)
  Q <- length(v)
  P <- (-1 + sqrt(1 + 8 * Q)) / 2
  M <- matrix(0,P,P)
  M[upper.tri(M,diag=T)] <- v
  M <- M + t(M)
  diag(M) <- diag(M) / 2
return(M)
}
#par.to.matrix(1:6)
