#' Helper functions for the penalized EM algorithm
#'
#' \code{vecInvDiag} is a helper function for quickly computing \eqn{(I + x \otimes y)^{-1}},
#' \code{tracePInvDiag} for quickly computing column sums of \eqn{(I + x \otimes y)^{-1}}. Both are used
#' in the penalized EM algorithm.
#'
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return for \code{vecInvDiag} a matrix defined by \eqn{(I + x \otimes y)^{-1}}, for
#' \code{tracePInvDiag} a matrix containing the column sums of \eqn{(I + x \otimes y)^{-1}}.
#'
#' @keywords internal
vecInvDiag <- function(x, y) {
  z <- sapply(X = x, FUN = function(x_i) {
    1 / (1 + x_i * y)
  })
  return(z)
}

#' @rdname vecInvDiag
tracePInvDiag <- function(x, y) {
  z <- sapply(X = x, FUN = function(x_i) {
    sum(1 / (1 + x_i * y))
  })
  return(z)
}
