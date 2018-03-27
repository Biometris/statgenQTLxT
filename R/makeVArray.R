#' Create array of (inverted) variance matrices
#'
#' As in Zhou and Stephens, supplement page 13, create the array of (inverted)
#' variance matrices for the transformed genotypes per individual as given by
#' equation 4.
#'
#' @param Vg a p x p symmetric matrix of genetic variance components
#' @param Ve a p x p symmetric matrix of environmental variance components
#' @param Dk a vector of length n containing the eigenvalues obtained by the
#' eigen-decomposition of the kinship matrix K: \eqn{K = Uk * Dk * t(Uk)}
#'
#' @return A list of n p x p matrices \eqn{v_l} where \eqn{v_l =
#' Dk_{l,l} * Vg + Ve \forall l = 1, ..., n}.\cr
#' When using \code{makeVInvArray} the output matrices are inverted.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409.
#'
#' @keywords internal
makeVArray <- function(Vg,
                       Ve,
                       Dk) {
  n <- length(Dk)
  p <- ncol(Vg)
  VArray <- array(dim = c(n, p, p))
  for (i in 1:n) {
    VArray[i, , ] <- as.matrix(Dk[i] * Vg + Ve)
  }
  return(VArray)
}

#' @rdname makeVArray
makeVInvArray <- function(Vg, Ve, Dk) {
  n <- length(Dk)
  p <- ncol(Vg)
  VInvArray <- array(dim = c(n, p, p))
  for (i in 1:n) {
    VInvArray[i, , ] <- solve(Dk[i] * Vg + Ve)
  }
  return(VInvArray)
}
