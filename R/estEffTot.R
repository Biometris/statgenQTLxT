## Helper function for estimating effect.
## Function is a wrapper around estEffsCPP, useful for usage in chromosome
## specific calculations.
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
  segMarkers <- which(abs(allFreq) >= 1 - MAF)
  ## Add snpCovariates to segregating markers.
  excludedMarkers <- union(c(segMarkers, ncol(markers) + 1),
                           exclMarkers(snpCov = snpCov, markers = markers,
                                       allFreq = allFreq))
  if (length(snpCov) > 0) {
    effEstSnpCovs <- lapply(X = snpCov, FUN = function(snp) {
      estEffsCPP(y = Y, w = X[, -which(colnames(X) == snp)],
                 x = as.matrix(markers[, snp, drop = FALSE]),
                 vg = VgRed, ve = VeRed, k = as.matrix(K),
                 estCom = estCom, nCores = nCores)
    })
    effEstSnpCov <- list(effs = do.call(cbind, lapply(X = effEstSnpCovs, `[[`, "effs")),
                         effsSe = do.call(cbind, lapply(X = effEstSnpCovs, `[[`, "effsSe")),
                         pVals = sapply(X = effEstSnpCovs, `[[`, "pVals"),
                         effsCom = sapply(X = effEstSnpCovs, `[[`, "effsCom"),
                         effsComSe = sapply(X = effEstSnpCovs, `[[`, "effsComSe"),
                         pValsCom = sapply(X = effEstSnpCovs, `[[`, "pValsCom"),
                         pValsQtlE = sapply(X = effEstSnpCovs, `[[`, "pValsQtlE"))
  } else {
    ## Set to NULL so binding can be done in next step.
    effEstSnpCov <- NULL
  }
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(markers)[-excludedMarkers]
  trtNames <- colnames(Y)
  effEst <- estEffsCPP(y = Y, w = X,
                       x = as.matrix(markers[, -excludedMarkers]),
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

