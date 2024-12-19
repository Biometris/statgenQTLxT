#' Functions for calculating kinship matrices
#'
#' A collection of functions for calculating kinship matrices using different
#' algorithms. The following algorithms are included: astle (Astle and Balding,
#' 2009), Identity By State (IBS) and VanRaden (VanRaden, 2008) for
#' marker matrices and multiAllKin for three-dimensional marker arrays as
#' created by \code{\link{readIBDProbs}}.
#'
#' @section Marker matrices:
#' In all algorithms the input matrix \code{X} is first cleaned, i.e. markers
#' with a variance of 0 are excluded from the calculation of the kinship matrix.
#' Then some form of scaling is done which differs per algorithm. This gives a
#' scaled matrix \code{Z}. The matrix \eqn{ZZ^t / denominator} is returned.
#' By default the denominator is equal to the number of columns in \code{Z} for
#' \code{astle} and \code{IBS} and \eqn{2 * p * (1-p)} where
#' \eqn{p = colSums(X) / (2 * nrow(X))} for \code{vanRaden}. This denominator
#' can be overwritten by the user, e.g. when computing kinship matrices by
#' splitting \code{X} in smaller matrices and then adding the results together
#' in the end.
#'
#' @section Three-dimensional marker arrays:
#' The kinship matrix K is computed as \eqn{\sum(Xm * Xm^t * dm) / \sum(dm)}
#' where Xm is the genotype x founders matrix for marker m and dm the sum
#' of the distances of marker m to its neighbours halved. Diagonals of K are set
#' to 1. Dividing by \eqn{\sum(dm)} can be overwritten by providing a value for
#' \code{denominator}.
#'
#' @param X An n x m marker matrix with genotypes in the rows (n) and markers in
#' the columns (m).
#' @param map An optional marker map. Only used when \code{X} is a 3 dimensional array.
#' @param method The method used for computing the kinship matrix. If \code{X} is a
#' matrix method cannot be "multiAllKin", if it is a three-dimensional array it
#' is automatically set to "multiAllKin".
#' @param denominator A numerical value. See details.
#'
#' @returns An n x n kinship matrix.
#'
#' @references Astle W., Balding D. J. (2009) Population structure and cryptic
#' relatedness in genetic association studies, Stat. Sci., November 2009,
#' Vol. 24, no. 4, p. 451–471.
#' @references VanRaden P.M. (2008) Efficient methods to compute genomic
#' predictions. J Dairy Sci, November 2008, Vol. 91 p. 4414–4423.
#'
#' @examples X <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1), nrow = 3)
#' kinship(X, method = "astle")
#' kinship(X, method = "IBS")
#' kinship(X, method = "vanRaden")
#'
#' @export
kinship <- function(X,
                    map = NULL,
                    method = c("astle", "IBS", "vanRaden", "multiAllKin"),
                    denominator = NULL) {
  method = match.arg(method);
  if (!is.null(denominator)) {
    chkNum(denominator, min = 0)
  }
  if (length(dim(X)) == 2 && method == "multiAllKin") {
    stop("method multiAllKin is only possible for three-dimensional array X.\n")
  } else if (length(dim(X)) == 3) {
    method <- "multiAllKin"
  }
  if (method == "multiAllKin") {
    ## Create an array of values for correction for position on the genome.
    ## Has to be done per chromosome since pos isn't necessary cumulative.
    posCor <- unlist(sapply(X = unique(map$chr), FUN = function(chr) {
      pos <- map[map$chr == chr, "pos"]
      ## First and last marker need special treatment. For those just take
      ## double the distance to the next/previous marker.
      posCor <- c(pos[2] - pos[1], diff(pos, lag = 2) / 2,
                  pos[length(pos)] - pos[length(pos) - 1])
    }))
    K <- multiAllKin(X, posCor, denominator)
    rownames(K) <- colnames(K) <- rownames(X)
  } else {
    K <- statgenGWAS::kinship(X = X, method = method, denominator = denominator)
  }
  return(K)
}
