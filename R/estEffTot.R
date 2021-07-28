## Helper function for estimating effect.
## Function is a wrapper around estEffsCPP, useful for usage in chromosome
## specific calculations.
#' @keywords internal
estEffTot <- function(markers,
                      map,
                      minCofactorProximity = 50,
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
  minDist <- minCofactorProximity / 2
  # segMarkers <- which(abs(allFreq) >= 1 - MAF)markersRed
  segMarkers <- which(abs(allFreq) < MAF)
  excludedMarkers <- c(segMarkers, ncol(markers) + 1)
  map <- map[-excludedMarkers, ]
  markers <- markers[, -excludedMarkers]
  if (length(snpCov) > 0) {
    ## Order snp cofactors by position on chromosome.
    snpCov <- rownames(map[rownames(map) %in% snpCov, ])
    snpCovPos <- map[snpCov, "pos"]
    snpRanges <- sapply(snpCovPos, FUN = function(pos) {
      map[["pos"]] > pos - minDist & map[["pos"]] < pos + minDist
    })
    effEstSnps <- lapply(X = colnames(markers), FUN = function(snp) {
      snpPos <- which(colnames(markers) == snp)
      snpCovSnp <- snpCov[snpRanges[snpPos, ]]
      effEstSnp <- estEffsCPP(y0 = Y,
                              w0 = X[, !colnames(X) %in% snpCovSnp, drop = FALSE],
                              x0 = as.matrix(markers[, snp, drop = FALSE]),
                              vg = Vg, ve = Ve, k = K, estCom = estCom,
                              nCores = nCores)
      effEstSnp$snp <- snp
      return(effEstSnp)
    })
    effEst <-
      list(effs = do.call(cbind, lapply(X = effEstSnps, `[[`, "effs")),
           effsSe = do.call(cbind, lapply(X = effEstSnps, `[[`, "effsSe")),
           pVals = unlist(sapply(X = effEstSnps, `[[`, "pVals")),
           effsCom = unlist(sapply(X = effEstSnps, `[[`, "effsCom")),
           effsComSe = unlist(sapply(X = effEstSnps, `[[`, "effsComSe")),
           pValsCom = unlist(sapply(X = effEstSnps, `[[`, "pValsCom")),
           pValsQtlE = unlist(sapply(X = effEstSnps, `[[`, "pValsQtlE")),
           snps = unlist(sapply(X = effEstSnps, `[[`, "snps")))
  } else {
    ## Compute all effects in one go.
    effEst <- estEffsCPP(y0 = Y, w0 = X,
                         x0 = as.matrix(markers), vg = Vg, ve = Ve, k = K,
                         estCom = estCom, nCores = nCores)
  }
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(markers)[-excludedMarkers]
  trtNames <- colnames(Y)
  ## Bind results together.
  pValues <- effEst$pVals
  effs <- effEst$effs
  effsSe <- effEst$effsSe
  pValCom <- effEst$pValsCom
  effsCom <- effEst$effsCom
  effsComSe <- effEst$effsComSe
  pValQtlE <- effEst$pValsQtlE
  names(pValues) <- colnames(effs) <- colnames(effsSe) <- names(pValCom) <-
    names(effsCom) <- names(effsComSe) <- names(pValQtlE) <- snpNames
  rownames(effs) <- rownames(effsSe) <- trtNames
  return(list(pValues = pValues, effs = effs, effsSe = effsSe,
              pValCom = pValCom, effsCom = effsCom, effsComSe = effsComSe,
              pValQtlE = pValQtlE))
}

