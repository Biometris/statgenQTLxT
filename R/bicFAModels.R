bicFAModels <- function(Y,
  K,
  X,
  faModels,
  maxDiag = 100,
  maxIter = 20000,
  tolerance = 1e-4,
  stopIfDecreasing = FALSE) {

  p <- ncol(Y)
  n <- nrow(Y)
  nc <- ncol(X)

  comparison <- sapply(X = 1:nrow(faModels), FUN = function(m) {
    EMFAResult <- EMFA(Y = Y,
      K = K,
      X = X,
      maxIter = maxIter,
      tolerance = tolerance,
      mG = faModels$mG[m],
      mE = faModels$mE[m],
      maxDiag = maxDiag,
      CmHet = faModels$hetG[m],
      DmHet = faModels$hetE[m],
      computeLogLik = TRUE,
      stopIfDecreasing = stopIfDecreasing)

    nPar <- (sum(faModels[m, ]) + nc) * p
    bic <- -2 * EMFAResult$logLik2 + log(n * p) * nPar
    aic <- -2 * EMFAResult$logLik2 + 2 * (n * p / (n * p - nPar - 1)) * nPar

    return(list(mG = faModels$m.G[m],
      mE = faModels$m.E[m],
      hetG = faModels$het.G[m],
      hetE = faModels$het.E[m],
      LL = EMFAResult$logLik2,
      LLE = EMFAResult$logLik,
      bic = bic,
      aic = aic,
      nIter = EMFAResult$n.iter,
      converged = EMFAResult$converged,
      decreased = EMFAResult$decreased))
  })
  return(as.data.frame(t(comparison)))
}
