#' Helper function for accessing parallel computing functions.
#'
#' @noRd
#' @keywords interanal
getOper <- function(x) {
  if (x) {
    `%dopar%`
  } else {
    `%do%`
  }
}
