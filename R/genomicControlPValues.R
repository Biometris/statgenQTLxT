#' Correction of p-values based on genomic inflation
#'
#' Correction of p-values based on the genomic inflation factor, as in Devlin and Roeder (1999). It is
#' assumed that the p-values come from an F-test with df1 = 1 and df2 = nObs - nCov - 2.
#'
#' @param pVals a numeric vector of p-values between 0 and 1; may contain NA's.
#' @param nObs an integer > 0 indicating the number of individuals.
#' @param nCov a positive integer indicating the number of covariables.
#' @return a list with two components:
#' \itemize{
#' \item{\code{pValues} a vector of p-values corrected by the genomic inflation factor, with the same NA's
#' as the input}.
#' \item{\code{inflation} the inflation factor}.
#' }
#'
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association studies. Biometrics,
#' December 1999, Vol. 55(4), p. 997-1004.

genomicControlPValues <- function(pVals,
  nObs,
  nCov = 0) {
  ## Check input.
  if (missing(pVals) || !is.numeric(pVals) || any(pVals, na.rm = TRUE) < 0 || any(pVals, na.rm = TRUE) > 1)
    stop("pVals should be a numeric vector with values between 0 and 1.\n")
  if (missing(nObs) || length(nObs) > 1 || !is.numeric(nObs) || nObs != round(nObs) || nObs < 1)
    stop("nObs should be a single integer > 0.\n")
  if (length(nCov) > 1 || !is.numeric(nCov) || nCov != round(nCov) || nCov < 0)
    stop("nCov should be a single positive integer.\n")

  ## Compute degree of freedom.
  df2 <- nObs - nCov - 2
  pValsNew <- pVals
  ## Compute F-values from input p-values.
  FVals <- qf(na.omit(pVals), df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute inflation factor as in Devlin and Roeder.
  inflation <- median(FVals, na.rm = TRUE) / qf(0.5, df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute new F-values and p-values.
  FValsNew <- FVals / inflation
  pValsNew[!is.na(pVals)] <- pf(FValsNew, df1 = 1, df2 = df2, lower.tail = FALSE)
  return(list(pValues = pValsNew, inflation = inflation))
}
