#' tracePInvDiag
#'
#' Helper function for quickly computing column sums of \eqn{(I + x \otimes y)^{-1}}. Used in the
#' penalized EM algorith.
#'
#' @inheritParams vecInvDiag
#'
#' @return a matrix containing the column sums of \eqn{(I + x \otimes y)^{-1}}
#'
#' @export
#'
#' @examples tracePInvDiag(1:2, 1:3)

tracePInvDiag <- function(x, y) {
  z <- sapply(x, function(x_i) {sum(1 / (1 + x_i * y))})
  return(z)
}
