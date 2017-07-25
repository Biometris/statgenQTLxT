#' kinshipTransform
#'
#' Helper function for computing ... Used for calculating kinship matrices.
#'
#' @param K an n by n matrix
#'
#' @return the ...

kinshipTransform <- function(K) {
  n <- ncol(K)
  return((sum(diag(K)) - rep(1, times = n) %*% K %*% rep(1, times = n) / n) / (n - 1))
}
