#' Helper function for estimating effect
#'
#' Function is a wrapper around estEffsCPP, useful for usage in chromosome
#' specific calculations.
#'
#' @keywords internal
estEffTot <- function(markers,
                      X,
                      Y,
                      K,
                      XRed,
                      Vg,
                      Ve,
                      VgRed,
                      VeRed,
                      snpCov,
                      allFreq,
                      MAF,
                      estCom,
                      nCores = NULL) {
  ## Determine segregating markers.
  segMarkers <- which(allFreq <= MAF | allFreq >= (1 - MAF))
  ## Add snpCovariates to segregating markers.
  excludedMarkers <- union(c(segMarkers, ncol(markers) + 1),
                           exclMarkers(snpCov = snpCov, markers = markers,
                                       allFreq = allFreq))
  if (length(snpCov) > 0) {
    effEstSnpCov <- estEffsCPP(y0 = Y, w0 = XRed,
                               x0 = as.matrix(markers[, snpCov, drop = FALSE]),
                               vg = VgRed, ve = VeRed, k = as.matrix(K),
                               estCom = estCom, nCores = nCores)
  } else {
    ## Set to NULL so binding can be done in next step.
    effEstSnpCov <- NULL
  }
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(markers)[-excludedMarkers]
  trtNames <- colnames(Y)
  effEst <- estEffsCPP(y0 = Y, w0 = X,
                       x0 = as.matrix(markers[, -excludedMarkers]),
                       vg = Vg, ve = Ve, k = as.matrix(K), estCom = estCom,
                       nCores = nCores)
  pValues <- c(effEst$pVals, effEstSnpCov$pVals)
  effs <- cbind(effEst$effs, effEstSnpCov$effs)
  effsSe <- cbind(effEst$effsSe, effEstSnpCov$effsSe)
  pValCom <- c(effEst$pValsCom, effEstSnpCov$pValsCom)
  effsCom <- c(effEst$effsCom, effEstSnpCov$effsCom)
  effsComSe <- c(effEst$effsComSe, effEstSnpCov$effsComSe)
  pValQtlE <- c(effEst$pValsQtlE, effEstSnpCov$pValsQtlE)
  names(pValues) <- colnames(effs) <- colnames(effsSe) <- names(pValCom) <-
    names(effsCom) <- names(effsComSe) <- names(pValQtlE) <-
    c(snpNames, snpCov)
  rownames(effs) <- rownames(effsSe) <- trtNames
  return(list(pValues = pValues, effs = effs, effsSe = effsSe,
              pValCom = pValCom, effsCom = effsCom, effsComSe = effsComSe,
              pValQtlE = pValQtlE))
}

