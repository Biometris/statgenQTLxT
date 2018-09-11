## Helper function for estimating effect.
## Function is a wrapper around estEffsCPP, usefull for usage in chromosome
## specific calculations.
#' @keywords internal
estEffTot <- function(markers,
                      X,
                      Y,
                      K,
                      XRed,
                      Vg,
                      Ve,
                      snpCov,
                      allFreq,
                      MAF,
                      estCom) {
  p <- ncol(Y)
  segMarkers <- which(allFreq < MAF | allFreq > 1 - MAF)
  ## Add snpCovariates to segregating markers.
  excludedMarkers <- union(c(segMarkers, ncol(markers) + 1),
                           exclMarkers(snpCov = snpCov, markers = markers,
                                       allFreq = allFreq))
  if (!is.null(snpCov)) {
    effEstSnpCov <- estEffsCPP(y = as.matrix(Y), w = as.matrix(XRed),
                               x = as.matrix(markers[, snpCov, drop = FALSE]),
                               vg = as.matrix(Vg), ve = as.matrix(Ve),
                               k = as.matrix(K), estCom = estCom)
  } else {
    ## Set to NULL so binding can be done in next step.
    effEstSnpCov <- NULL
  }
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(markers)[-excludedMarkers]
  trtNames <- colnames(Y)
  effEst <- estEffsCPP(y = as.matrix(Y), w = as.matrix(X),
                       x = as.matrix(markers[, -excludedMarkers]),
                       vg = as.matrix(Vg), ve = as.matrix(Ve),
                       k = as.matrix(K), estCom = estCom)
  FVals <- c(effEst$fVals, effEstSnpCov$fVals)
  pValues <- pf(q = FVals, df1 = p, df2 = effEst$dfFull, lower.tail = FALSE)
  effs <- cbind(effEst$effs, effEstSnpCov$effs)
  effsSe <- cbind(effEst$effsSe, effEstSnpCov$effsSe)
  FValCom <- c(effEst$fValCom, effEstSnpCov$fValCom)
  pValCom <- pf(q = FValCom, df1 = 1, df2 = effEst$dfCom, lower.tail = FALSE)
  effsCom <- c(effEst$effsCom, effEstSnpCov$effsCom)
  effsComSe <- c(effEst$effsComSe, effEstSnpCov$effsComSe)
  FValQtlE <- c(effEst$fValQtlE, effEstSnpCov$fValQtlE)
  pValQtlE <- pf(q = FValQtlE, df1 = p - 1, df2 = effEst$dfFull,
                 lower.tail = FALSE)
  names(pValues) <- colnames(effs) <- colnames(effsSe) <- names(pValCom) <-
    names(effsCom) <- names(effsComSe) <- names(pValQtlE) <-
    c(snpNames, snpCov)
  rownames(effs) <- rownames(effsSe) <- trtNames
  return(list(pValues = pValues, effs = effs, effsSe = effsSe,
              pValCom = pValCom, effsCom = effsCom, effsComSe = effsComSe,
              pValQtlE = pValQtlE))
}

