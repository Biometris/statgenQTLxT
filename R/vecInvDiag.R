#' vecInvDiag
#'
#' Helper function for quickly computing \eqn{(I + x \otimes y)^{-1}}. Used in the penalized EM algorith.
#'
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return a matrix defined by \eqn{(I + x \otimes y)^{-1}}
#'
#' @examples vecInvDiag(1:2, 1:3)

vecInvDiag <- function(x, y) {
  z <- sapply(x, function(x_i) {1 / (1 + x_i * y)})
  return(z)
}
