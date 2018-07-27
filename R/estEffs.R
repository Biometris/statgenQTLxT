#' Estimates for covariates
#'
#' Compute the estimates and standard errors for the covariates in the input
#' matrix W.
#'
#' @param Y An n x p matrix of observed phenotypes, on p traits or environments
#' for n genotypes. No missing values are allowed.
#' @param W An n x c covariate matrix, c being the number of covariates and n
#' being the number of genotypes. c has to be at least one (typically an
#' intercept). No missing values are allowed.
#' @param X An n x ns matrix of marker scores. Neither missing values nor
#' non-segregating markers are allowed.
#' @param Vg A p x p matrix of genetic covariances.
#' @param Ve A p x p matrix of environmental covariances.
#' @param K An n x n genetic relatedness matrix.
#' @param returnSe Should standard errors and p-values be returned?
#' @param estCom Should the common SNP-effect model be fitted?
#' @param nChunks An integer, the number of blocks in which the calculations
#' should be split.
#'
#' @return A list containing the estimates, optionally the standard errors of
#' the estimates and corresponding p-values. If \code{estCom = TRUE} also
#' common SNP-effects, their standard errors and corresponding p-values and
#' the p-values for QtlxE are output.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal
estEffs <- function(Y,
                    W = matrix(data = 1, nrow = nrow(Y)),
                    X,
                    Vg,
                    Ve,
                    K,
                    returnSe = TRUE,
                    estCom = FALSE,
                    nChunks = min(ncol(X), ceiling(ncol(Y) * ncol(W) *
                                                     ncol(X) / 50000))) {
  ## Y, W and X might be from the Matrix class. This function needs standard
  ## matrices for its computations.
  Y <- as.matrix(Y)
  W <- as.matrix(W)
  X <- as.matrix(X)
  ## Extract number of traits, individuals, covariates and SNPs.
  p <- ncol(Y)
  n <- nrow(Y)
  nc <- ncol(W)
  ns <- ncol(X)
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(X)
  genoNames <- colnames(Y)
  ## Compute eigen-decomposition of K.
  w <- eigen(K, symmetric = TRUE)
  Dk <- w$values
  Dk[Dk < max(Dk) * 1e-5] <- max(Dk) * 1e-5
  Uk <- w$vectors
  ## Transform Y, W and X.
  Y <- t(Y) %*% Uk
  W <- t(W) %*% Uk
  X <- t(X) %*% Uk
  ## Square each element of X.
  X2 <- X ^ 2
  # Define vInvArr
  vInvArr <- makeVInvArray(Vg = Vg, Ve = Ve , Dk = Dk)
  ## create a second instance of vInvArr, in matrix form.
  vInvArrRed <- vInvArr
  dim(vInvArrRed) <- c(n, p ^ 2)
  ## Compute quantities that are independent of the SNPs.
  VBeta <- matrix(rowSums(sapply(X = 1:n, FUN = function(i) {
    kronecker(tcrossprod(W[, i]), vInvArr[i, , ])
  })), ncol = p * nc)
  v <- rowSums(sapply(X = 1:n, FUN = function(i) {
    kronecker(W[, i], vInvArr[i, , ] %*% Y[, i])
  }))
  ## VInvY is used is several equations so is computed once here.
  VInvY <- sapply(X = 1:n, FUN = function(i) {
    vInvArr[i, , ] %*% Y[, i]
  })
  ## Define output for effects.
  Eff <- matrix(data = 0, nrow = p, ncol = ns,
                dimnames = list(genoNames, snpNames))
  if (returnSe) {
    ## Define output for standard error of effects.
    EffSe <- Eff
    ## Compute SS0 for null model with the trait specific means only.
    est0 <- solve(VBeta, v)
    fitMean0 <- matrix(est0, ncol = length(est0) / p) %*% W
    SS0 <- LLQuadFormDiag(Y = Y - fitMean0, vInvArr = vInvArr)
    ## Compute VQ
    VQ <- numeric(p * nc)
    for (c in 1:nc) {
      VQ[(c - 1) * p + 1:p] <- VInvY %*% W[c, ]
    }
    ## Compute scalar part of SS1.
    qScal <- sum(Y * VInvY)
    ## Define output for SS1.
    SS1 <- setNames(rep(x = NA, times = ns), snpNames)
  }
  if (estCom) {
    ## Define output for common effects.
    EffCom <- setNames(numeric(ns), snpNames)
    ## Compute chunk independent quantities.
    s1 <- colSums(VInvY)
    s2 <- apply(X = vInvArr, MARGIN = 1, FUN = sum)
    S3 <- t(apply(X = vInvArr, MARGIN = 1, FUN = rowSums))
    if (returnSe) {
      ## Define outputs for SE of common effects and SS for common effects.
      EffSeCom <- EffCom
      SS1Com <- SS1
    }
  }
  ## Calculate chunk sizes.
  chunks <- split(1:ns, c(rep(1:nChunks, each = ns %/% nChunks),
                          rep(nChunks, each = ns %% nChunks)))
  ## All calculations are done for chunks of markers to avoid memory problems.
  for (ch in chunks) {
    ## Compute chunk length.
    nsCh <- length(ch)
    ## Compute position of first SNP in current chunk -used for filling results.
    snpPos0 <- ch[1] - 1
    ## Define SNP-dependent quantities.
    X2VinvX1Arr <- array(data = 0, dim = c(nc, nsCh, p ^ 2))
    ## Fill X2VinvX1Arr by looping over covariates.
    for (cv in 1:nc) {
      X2VinvX1Arr[cv, , ] <- X[ch, ] %*% (vInvArrRed * as.numeric(W[cv, ]))
    }
    ## Compute VBeta and V per SNP - computation only depends on individuals.
    VBetaSnp <- X2[ch, ] %*% vInvArrRed
    VSnp <- tcrossprod(X[ch, ], VInvY)
    if (returnSe) {
      ## Compute VSnpQ by looping over individuals.
      VSnpQ <- matrix(data = 0, nrow = nsCh, ncol = p)
      for (i in 1:n) {
        VSnpQ <- VSnpQ + X[ch, i] %o% VInvY[, i]
      }
    }
    ## The remaining calculations are SNP-dependent.
    for (snp in 1:nsCh) {
      ## Compute inverse of VBetaSnp.
      VBetaSnpInv <- solve(matrix(data = VBetaSnp[snp, ], ncol = p))
      ## Extract X2VinvX1 from array.
      X2VinvX1 <- matrix(t(X2VinvX1Arr[, snp, ]), ncol = p, byrow = TRUE)
      ## Compute inverse of XS.
      XSInv <- solve(VBeta - tcrossprod(X2VinvX1 %*% VBetaSnpInv, X2VinvX1))
      ## Compute SNP dependent effect.
      EffCovSnp <- XSInv %*% (v - X2VinvX1 %*% VBetaSnpInv %*% VSnp[snp, ])
      ## Compute SNP effect.
      Eff[, snpPos0 + snp] <- VBetaSnpInv %*%
        (VSnp[snp, ] - crossprod(X2VinvX1, EffCovSnp))
      if (returnSe) {
        ## Compute inverse of QSnp.
        QSnpInv <- solve(rbind(cbind(VBeta, X2VinvX1),
                               cbind(t(X2VinvX1), solve(VBetaSnpInv))))
        ## Compute SE of SNP effect.
        EffSe[, snpPos0 + snp] <- sqrt(diag(QSnpInv)[-(1:(p * nc))])
        ## Compute SS1 per SNP.
        QSnp <- c(VQ, VSnpQ[snp, ])
        SS1[snpPos0 + snp] <- qScal - crossprod(QSnp, QSnpInv %*% QSnp)
      } # End returnSe
    } # End loop over SNPs
    if (estCom) {
      ## Calculations for common effect are similar to those above.
      ## However the dimensions are lower and therefore some calculations
      ## have been simplified.
      ## Define SNP-dependent quantities.
      X2VinvX1ArrCom <- array(data = 0, dim = c(nc, nsCh, p))
      ## Fill X2VinvX1ArrCom by looping over covariates.
      for (cv in 1:nc) {
        X2VinvX1ArrCom[cv, , ] <- X[ch, ] %*% (S3 * as.numeric(W[cv, ]))
      }
      ## Compute common VBeta and V per SNP - only depends on individuals.
      VSnpCom <- X[ch, ] %*% s1
      VBetaSnpCom <- X2[ch, ] %*% s2
      if (returnSe) {
        ## Compute VSnpQ for common effect.
        vSnpQCom <- rowSums(X[ch, ] %*% colSums(VInvY))
      }
      ## The remaining calculations are SNP-dependent.
      for (snp in 1:nsCh) {
        ## Compute VBetaSnpInvCom - scalar value
        VBetaSnpInvCom <- 1 / VBetaSnpCom[snp]
        ## Extract X2VinvX1Com from array.
        X2VinvX1Com <- as.numeric(t(X2VinvX1ArrCom[, snp, ]))
        ## Compute inverse of XS for common effect.
        XSInvCom <- solve(VBeta - VBetaSnpInvCom * tcrossprod(X2VinvX1Com))
        ## Compute SNP dependent common effect.
        EffCovComSnp <- XSInvCom %*%
          (v - VBetaSnpInvCom * VSnpCom[snp] * X2VinvX1Com)
        ## Compute common effect.
        EffCom[snpPos0 + snp] <- VBetaSnpInvCom *
          (VSnpCom[snp] - X2VinvX1Com %*% EffCovComSnp)
        if (returnSe) {
          ## Compute inverse of QSnp for common effect.
          QSnpInvCom <- solve(rbind(cbind(VBeta, X2VinvX1Com),
                                    c(X2VinvX1Com, 1 / VBetaSnpInvCom)))
          ## Compute SE of common SNP effect.
          EffSeCom[snpPos0 + snp] <- sqrt(diag(QSnpInvCom)[-(1:(p * nc))])
          ## Compute SS1 for common effect.
          QSnpCom <- c(VQ, vSnpQCom[snp])
          SS1Com[snpPos0 + snp] <- qScal -
            crossprod(QSnpCom, QSnpInvCom %*% QSnpCom)
        } # End returnSe
      } # End loop over SNPs
    } # End estCom
  } # End loop over chuncks
  if (returnSe) {
    ## Compute degrees of freedom for full model.
    dfFull <- (n - nc - 1) * p
    ## Compute F-values and p-values for SNP effect.
    FVals <- ((SS0 - SS1) / SS1) * dfFull / p
    pVals <- pf(q = FVals, df1 = p, df2 = dfFull, lower.tail = FALSE)
    if (estCom) {
      ## Compute degrees of freedom for common effect model.
      dfCom  <- (n - nc) * p - 1
      ## Compute F-values and p-values for common SNP effect.
      FValsCom <- (SS0 - SS1Com) / SS1Com * dfCom
      pValsCom <- pf(q = FValsCom, df1 = 1, df2 = dfCom, lower.tail = FALSE)
      ## Compute F-values and p-values for QTLxE effect.
      FValsQtlE <- ((SS1Com - SS1) / SS1) * dfFull / (p - 1)
      pValsQtlE <- pf(q = FValsQtlE, df1 = p - 1, df2 = dfFull,
                      lower.tail = FALSE)
    }
  }
  ## Construct output.
  out <- list(effs = Eff)
  if (returnSe) {
    out$effsSe <- EffSe
    out$pVals <- pVals
  }
  if (estCom) {
    out$effsCom <- EffCom
    if (returnSe) {
      out$effsComSe <- EffSeCom
      out$pValCom <- pValsCom
      out$pValQtlE <- pValsQtlE
    }
  }
  return(out)
}

## Helper function for estimating effect.
## Function is a wrapper around estEffs, usefull for usage in chromosome
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
  segMarkers <- which(allFreq < MAF | allFreq > 1 - MAF)
  ## Add snpCovariates to segregating markers.
  excludedMarkers <- union(c(segMarkers, ncol(markers) + 1),
                           exclMarkers(snpCov = snpCov, markers = markers,
                                       allFreq = allFreq))
  if (!is.null(snpCov)) {
    effEstSnpCov <- estEffs(Y = Y, W = XRed,
                            X = markers[, snpCov, drop = FALSE], Vg = Vg,
                            Ve = Ve, K = K, estCom = estCom)
  } else {
    ## Set to NULL so binding can be done in next step.
    effEstSnpCov <- NULL
  }
  effEst <- estEffs(Y = Y, W = X, X = markers[, -excludedMarkers], Vg = Vg,
                    Ve = Ve, K = K, estCom = estCom)
  pValues <- c(effEst$pVals, effEstSnpCov$pVals)
  effs <- cbind(effEst$effs, effEstSnpCov$effs)
  effsSe <- cbind(effEst$effsSe, effEstSnpCov$effsSe)
  pValCom <- c(effEst$pValCom, effEstSnpCov$pValCom)
  effsCom <- c(effEst$effsCom, effEstSnpCov$effsCom)
  effsComSe <- c(effEst$effsComSe, effEstSnpCov$effsComSe)
  pValQtlE <- c(effEst$pValQtlE, effEstSnpCov$pValQtlE)
  return(list(pValues = pValues, effs = effs, effsSe = effsSe,
              pValCom = pValCom, effsCom = effsCom, effsComSe = effsComSe,
              pValQtlE = pValQtlE))
}

