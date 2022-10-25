#' Row bind data.frames
#'
#' Helper function for row binding data.frames with diffent columns.
#'
#' @param dfList A list of data.frames.
#'
#' @noRd
#' @keywords internal
dfBind <- function(dfList) {
  ## Filter empty data.frames from dfList
  dfList <- Filter(f = function(x) nrow(x) > 0, x = dfList)
  if (length(dfList) == 0) {
    return(data.frame())
  }
  ## Get variable names from all data.frames.
  allNms <- unique(unlist(lapply(dfList, names)))
  ## rbind all data.frames setting values for missing columns to NA.
  do.call(rbind,
          c(lapply(X = dfList, FUN = function(x) {
            nwDat <- sapply(X = setdiff(allNms, names(x)), FUN = function(y) {
              NA
            })
            data.frame(c(x, nwDat), stringsAsFactors = FALSE)
          }), make.row.names = FALSE)
  )
}

## Helper function for accessing parallel computing functions.
getOper <- function(x) {
  if (x) {
    `%dopar%`
  } else {
    `%do%`
  }
}
