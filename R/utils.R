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

#' Read IBD probabilities
#'
#' Read a file with IBD probabilities computed by the RABBIT software package.
#'
#' @param infile A character string, a link to a .csv file with IBD
#' probabilities.
#'
#' @return A list with two components:
#' \itemize{
#' \item{markArr} {A three-dimensional array containing the IBD probabilities
#' with genotypes in the rows, markers in the columns and founders in the
#' third dimension.}
#' \item{map} {A data.frame with two columns, chr(omosome) and pos(ition) with
#' genotypes as rownames.}
#' }
#' Both markArr and map can be used directly in \code{\link{createGData}}.
#'
#' @export
readIBDProbs <- function(infile) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 || tools::file_ext(infile) != "csv") {
    stop("infile should be a character string indicating a readable .csv file")
  }
  ## Read map and marker probabilities.
  markMap <- data.table::fread(infile, skip = "haploprob", fill = TRUE,
                               data.table = FALSE)
  ## Extract map.
  map <- data.frame(chr = as.numeric(markMap[3, -1]),
                    pos = as.numeric(markMap[4, -1]),
                    row.names = as.character(markMap[2, -1]))
  ## Get names of genotypes and compute number of founder alleles per genotype.
  genoNames <- unique(sapply(X = strsplit(x = markMap[5:nrow(markMap), 1],
                                          split = "_haplotype"), FUN = "[[", 1))
  nAlleles = (nrow(markMap) - 4) / length(genoNames)
  ## Convert markers to 3D array.
  markArr <- array(dim = c(length(genoNames), nrow(map), nAlleles))
  for (i in 1:nrow(map)) {
    markArr[, i, ] <- matrix(as.numeric(markMap[5:nrow(markMap), i + 1]),
                             ncol = nAlleles, byrow = TRUE)
  }
  ## Read founder names from file.
  foundNames <- data.table::fread(infile, header = FALSE, nrows = nAlleles + 2,
                                  skip = "haplotypes in order",
                                  data.table = FALSE, select = 3)
  foundNames <- as.character(foundNames[3:nrow(foundNames), ])
  ## Add dimnames to markers: genotypes x markers x founders.
  dimnames(markArr) <- list(genoNames, rownames(map), foundNames)
  return(list(markArr = markArr, map = map))
}

