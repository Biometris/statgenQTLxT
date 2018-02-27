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
#' common SNP-effects and their standard errors are output.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal
estimateEffects <- function(Y,
                            W = matrix(data = 1, nrow = nrow(Y)),
                            X,
                            Vg,
                            Ve,
                            K,
                            returnSe = TRUE,
                            estCom = FALSE,
                            nChunks = 200) {
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
  VBeta <- matrix(rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(tcrossprod(W[, m]), vInvArr[m, , ])
  })), ncol = p * nc)
  v <- rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(W[, m], vInvArr[m, , ] %*% Y[, m])
  }))
  ## Calculate chunk sizes.
  chunks <- split(1:ns, c(rep(1:nChunks, each = ns %/% nChunks),
                          rep(nChunks, each = ns %% nChunks)))
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
    ## Define output for SS1.
    SS1 <- setNames(rep(x = NA, times = ns), snpNames)
  }
  if (estCom) {
    ## Define output for common effects.
    EffCom <- setNames(numeric(ns), snpNames)
    ## Compute chunk independent quantities.
    s1 <- sapply(1:n, function(i) {
      sum(vInvArr[i, , ] %*% Y[, i])
    })
    s2 <- apply(vInvArr, MARGIN = 1, sum)
    S3 <- apply(vInvArr, MARGIN = 1, rowSums)
  }
  if (returnSe) {
    EffSeCom <- EffCom
    SS1Com <- SS1
  }
  for (ch in chunks) {
    nsCh <- length(ch)
    snpPos0 <- ch[1] - 1
    ## Define SNP-dependent quantities.
    VSnp <- matrix(data = 0, nrow = nsCh, ncol = p)
    VBetaSnp <- matrix(data = 0, nrow = nsCh, ncol = p ^ 2)
    X2VinvX1Arr <- array(data = 0, dim = c(nc, nsCh, p ^ 2))
    if (returnSe) {
      EffCov <- matrix(data = 0, nrow = nc * p, ncol = nsCh)
      QInv <- array(dim = c(nsCh, (nc + 1) * p, (nc + 1) * p))
    }
    ## 'Fill' the SNP-dependent quantities, looping over the individuals.
    for (i in 1:n) {
      VSnp <- VSnp + tcrossprod(X[ch, i], vInvArr[i, , ] %*% Y[, i])
      VBetaSnp <- VBetaSnp + tcrossprod(X2[ch, i], vInvArrRed[i, ])
      for (cv in 1:nc) {
        X2VinvX1Arr[cv, , ] <- X2VinvX1Arr[cv, , ] +
          W[cv, i] * tcrossprod(X[ch, i], vInvArrRed[i, ])
      }
    }
    ## Now do the remaining calculations, looping over the SNPs.
    for (snp in 1:nsCh) {
      VBetaSnpInv <- solve(matrix(data = VBetaSnp[snp, ], ncol = p))
      X2VinvX1 <- matrix(t(X2VinvX1Arr[, snp, ]), ncol = p, byrow = TRUE)
      XSInv <- solve(VBeta - tcrossprod(X2VinvX1 %*% VBetaSnpInv, X2VinvX1))
      EffCovSnp <- XSInv %*% (v - X2VinvX1 %*% VBetaSnpInv %*% VSnp[snp, ])
      Eff[, snpPos0 + snp] <- VBetaSnpInv %*%
        (VSnp[snp, ] - crossprod(X2VinvX1, EffCovSnp))
      if (returnSe) {
        EffCov[, snp] <- EffCovSnp
        QSnpInv <- solve(rbind(cbind(VBeta, X2VinvX1),
                               cbind(t(X2VinvX1), solve(VBetaSnpInv))))
        EffSe[, snpPos0 + snp] <- sqrt(diag(QSnpInv)[-(1:(p * nc))])
        QInv[snp, , ] <- QSnpInv
      }
    }
    if (returnSe) {
      qScal <- numeric(nsCh)
      VQ <- matrix(data = 0, nrow = nsCh, ncol = p * nc)
      VSnpQ <- matrix(data = 0, nrow = nsCh, ncol = p)
      for (i in 1:n) {
        Res <- t(matrix(data = 1, nrow = p, ncol = nsCh) * Y[, i]) -
          X[ch, i] * t(Eff[, ch]) - W[1, i] * t(EffCov[1:p, ])
        if (nc > 1) {
          for (ec in 2:nc) {
            Res <- Res - W[ec, i] * t(EffCov[(ec - 1) * p + 1:p, ])
          }
        }
        ResVInv <- Res %*% vInvArr[i, , ]
        qScal <- qScal + rowSums(Res * ResVInv)
        VSnpQ <- VSnpQ + X[ch, i] * ResVInv
        VQ[, 1:p] <- VQ[, 1:p] + W[1, i] * ResVInv
        if (nc > 1) {
          for (ec in 2:nc) {
            VQ[, (ec - 1) * p + 1:p] <- VQ[, (ec - 1) * p + 1:p] +
              W[ec, i] * ResVInv
          }
        }
      }
      for (snp in 1:nsCh) {
        QSnp <- matrix(c(VQ[snp, ], VSnpQ[snp, ]))
        SS1[snpPos0 + snp] <- qScal[snp] -
          crossprod(QSnp, QInv[snp, , ] %*% QSnp)
      }
    }
    if (estCom) {
      ## Define SNP-dependent quantities
      X2VinvX1ArrCom <- array(data = 0, dim = c(nc, nsCh, p))
      if (returnSe) {
        EffCovCom <- matrix(data = 0, nrow = nc * p, ncol = nsCh)
        QInvCom <- array(dim = c(nsCh, nc * p + 1, nc * p + 1))
      }
      ## 'Fill' the SNP-dependent quantities, looping over the individuals.
      VSnpCom <- X[ch, ] %*% s1
      VBetaSnpCom <- X2[ch, ] %*% s2
      for (i in 1:n) {
        for (cv in 1:nc) {
          X2VinvX1ArrCom[cv, , ] <- X2VinvX1ArrCom[cv, , ] +
            W[cv, i] * tcrossprod(X[ch, i], S3[, i])
        }
      }
      for (snp in 1:nsCh) {
        VBetaSnpInvCom <- 1 / VBetaSnpCom[snp]
        X2VinvX1Com <- matrix(data = X2VinvX1ArrCom[, snp, ], ncol = p,
                              byrow = TRUE)
        XSInvCom <- solve(VBeta - VBetaSnpInvCom * crossprod(X2VinvX1Com))
        EffCovComSnp <- XSInvCom %*%
          (v - VBetaSnpInvCom * VSnpCom[snp] * t(X2VinvX1Com))
        EffCom[snpPos0 + snp] <- VBetaSnpInvCom %*%
          (VSnpCom[snp] - X2VinvX1Com %*% EffCovComSnp)
        if (returnSe) {
          EffCovCom[, snp] <- EffCovComSnp
          QSnpInvCom <- solve(rbind(cbind(VBeta, t(X2VinvX1Com)),
                                    c(X2VinvX1Com, 1 / VBetaSnpInvCom)))
          EffSeCom[snpPos0 + snp] <- sqrt(diag(QSnpInvCom)[-(1:(p * nc))])
          QInvCom[snp, , ] <- QSnpInvCom
        }
      }
      if (returnSe) {
        ### Null model with the trait specific means only
        #    (which should be the model for which the Vg and Ve estimates were obtained)
        qScalCom <- VSnpQCom <- numeric(nsCh)
        VQCom <- matrix(data = 0, nrow = nsCh, ncol = p * nc)
        for (i in 1:n) {
          ResCom <- t(matrix(data = 1, nrow = p, ncol = nsCh) * Y[, i]) -
            X[ch, i] %*% matrix(data = 1, ncol = p) * EffCom[ch] -
            W[1, i] * t(EffCovCom[1:p, ])
          if (nc > 1) {
            for (ec in 2:nc) {
              ResCom <- ResCom - t(EffCovCom[(ec - 1) * p + 1:p, ] * W[ec, i])
            }
          }
          ResVInvCom <- ResCom %*% vInvArr[i, , ]
          qScalCom <- qScalCom + rowSums(ResCom * ResVInvCom)
          VSnpQCom <- VSnpQCom + X[ch, i] * rowSums(ResVInvCom)
          VQCom[, 1:p] <- VQCom[, 1:p] + W[1, i] * ResVInvCom
          if (nc > 1) {
            for (ec in 2:nc) {
              VQCom[, (ec - 1) * p + 1:p] <- VQCom[, (ec - 1) * p + 1:p] +
                W[ec, i] * ResCom %*% vInvArr[i, , ]
            }
          }
        }
        for (snp in 1:nsCh) {
          QSnpCom <- matrix(c(VQCom[snp, ], VSnpQCom[snp]))
          SS1Com[snpPos0 + snp] <- qScalCom[snp] -
            crossprod(QSnpCom, QInvCom[snp, , ] %*% QSnpCom)
        }
      }
    } # End estCom
  } # End loop over chuncks
  if (returnSe) {
    ## Compute degrees of freedom
    dfFull <- (n - nc - 1) * p
    FVals <- ((SS0 - SS1) / SS1) * dfFull / p
    pVals <- pf(q = FVals, df1 = p, df2 = dfFull, lower.tail = FALSE)
    if (estCom) {
      dfCom  <- (n - nc) * p - 1
      FValCom <- (SS0 - SS1Com) / SS1Com * dfCom
      pValCom <- pf(q = FValCom, df1 = 1, df2 = dfCom, lower.tail = FALSE)
      FValQtlE <- ((SS1Com - SS1) / SS1) * dfFull / (p - 1)
      pValQtlE <- pf(q = FValQtlE, df1 = p - 1, df2 = dfFull,
                     lower.tail = FALSE)
    }
  }
  ## Construct output.
  out <- list(effects = Eff)
  if (returnSe) {
    out$effectsSe <- EffSe
    out$pVals <- pVals
  }
  if (estCom) {
    out$effectsCom <- EffCom
    if (returnSe) {
      out$effectsComSe <- EffSeCom
      out$pValCom <- pValCom
      out$PValQtlE <- pValQtlE
    }
  }
  return(out)
}
