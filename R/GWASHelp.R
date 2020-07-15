#' Select markers to be excluded from GWAS scan.
#'
#' Helper function for selecting markers to be excluded from GWAS scan.
#' Markers are excluded if they are identical to any of the snpCovariates
#' (including the snpCovariates themselves).
#'
#' @param snpCov A character vector of snpCovariates.
#' @param markers A matrix with marker information.
#' @param allFreq A numerical vector of allele frequencies of the markers in
#' \code{markers}. This could be computed from markers as well but it is
#' needed in the general algorithm so to not redo things unnecessarily it is
#' not redone here.
#'
#' @return A numerical vector of markers to be exluded from the GWAS scan.
#'
#' @keywords internal
exclMarkers <- function(snpCov,
                        markers,
                        allFreq,
                        ref = NULL) {
  exclude <- integer()
  if (any(snpCov %in% colnames(markers))) {
    snpCovNumbers <- which(colnames(markers) %in% snpCov)
    if (length(dim(markers)) == 2) {
      for (snp in snpCovNumbers) {
        ## Rough selection based on allele frequency. Done for speed.
        candidates <- which(allFreq == allFreq[snp])
        ## Exclude all snps that are identical to snps in snpCovariates.
        snpInfo <- as.numeric(markers[, snp])
        exclude <- union(exclude,
                         candidates[apply(X = markers[, candidates,
                                                      drop = FALSE],
                                          MARGIN = 2, FUN = function(x) {
                                            identical(as.numeric(x), snpInfo)
                                          })])
      }
    } else if (length(dim(markers)) == 3) {
      ## Compute mean value for reference allele.
      allMeans <- apply(markers[ , , -ref], c(3, 2), mean)
      for (snp in snpCovNumbers) {
        for (allele in rownames(allMeans[allMeans[, snp] != 0, ])) {
          ## Rough selection based on mean. Done for speed.
          candidates <- which(allMeans == allMeans[allele, snp], arr.ind = TRUE)
          exclude <- union(exclude,
                           candidates[apply(X = candidates, MARGIN = 1,
                                            FUN = function(m) {
                                              identical(markers[, m[2], m[1]],
                                                        markers[, snp, allele])
                                            }), 2])
        }
      }
      ## Rough selection based on mean. Done for speed.
      candidates <- which(allMeans == allMeans[snp])
      ## Exclude all snps that are identical to snps in snpCovariates.
      snpInfo <- as.numeric(markers[, snp, ])
      exclude <- union(exclude,
                       candidates[apply(X = markers[, candidates, ,
                                                    drop = FALSE],
                                        MARGIN = 2, FUN = function(x) {
                                          identical(as.numeric(x), snpInfo)
                                        })])

    }
  }
  return(exclude)
}

#' Correction of p-values based on genomic inflation
#'
#' Correction of p-values based on the genomic inflation factor, as in Devlin
#' and Roeder (1999). It is assumed that the p-values come from an F-test with
#' df1 = 1 and df2 = nObs - nCov - 2.
#'
#' @param pVals A numeric vector of p-values between 0 and 1; may contain NA's.
#' @param nObs An integer > 0 indicating the number of individuals.
#' @param nCov An integer > 0 indicating the number of covariables.
#'
#' @return A list with two components:
#' \itemize{
#' \item{\code{pValues} a vector of p-values corrected by the genomic inflation
#' factor, with the same NA's as the input}.
#' \item{\code{inflation} the inflation factor}.
#' }
#'
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association
#' studies. Biometrics, December 1999, Vol. 55(4), p. 997-1004.
#'
#' @keywords internal
genCtrlPVals <- function(pVals,
                         nObs,
                         nCov = 0) {
  ## Check input.
  if (missing(pVals) || !is.numeric(pVals) || any(pVals < 0, na.rm = TRUE) ||
      any(pVals > 1, na.rm = TRUE)) {
    stop("pVals should be a numeric vector with values between 0 and 1.\n")
  }
  if (missing(nObs) || length(nObs) > 1 || !is.numeric(nObs) ||
      nObs != round(nObs) || nObs < 1) {
    stop("nObs should be a single positive integer.\n")
  }
  if (length(nCov) > 1 || !is.numeric(nCov) || nCov != round(nCov) ||
      nCov < 0) {
    stop("nCov should be a single non negative integer.\n")
  }
  ## Compute degree of freedom.
  df2 <- nObs - nCov - 2
  pValsNew <- pVals
  ## Compute F-values from input p-values.
  fVals <- qf(p = na.omit(pVals), df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute inflation factor as in Devlin and Roeder.
  inflation <- median(fVals, na.rm = TRUE) /
    qf(p = 0.5, df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute new F-values and p-values.
  fValsNew <- fVals / inflation
  pValsNew[!is.na(pVals)] <- pf(q = fValsNew, df1 = 1, df2 = df2,
                                lower.tail = FALSE)
  return(list(pValues = pValsNew, inflation = inflation))
}
