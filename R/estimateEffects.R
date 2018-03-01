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
                            nChunks = 50) {
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
    ## Compute scalar part.
    qScal <- sum(sapply(X = 1:n, FUN = function(i) {
      as.numeric(crossprod(Y[, i], vInvArr[i, ,] %*% Y[, i]))
    }))
    VQ <- numeric(p * nc)
    YSnpVInv <- sapply(X = 1:n, FUN = function(i) {
      vInvArr[i, , ] %*% Y[, i]
    })
    for (i in 1:n) {
      for (c in 1:nc) {
        VQ[(c - 1) * p + 1:p] <- VQ[(c - 1) * p + 1:p] + W[c, i] * YSnpVInv[, i]
      }
    }
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
    if (returnSe) {
      EffSeCom <- EffCom
      SS1Com <- SS1
    }
  }
  ## Calculate chunk sizes.
  chunks <- split(1:ns, c(rep(1:nChunks, each = ns %/% nChunks),
                          rep(nChunks, each = ns %% nChunks)))
  for (ch in chunks) {
    nsCh <- length(ch)
    snpPos0 <- ch[1] - 1
    ## Define SNP-dependent quantities.
    VSnp <- matrix(data = 0, nrow = nsCh, ncol = p)
    X2VinvX1Arr <- array(data = 0, dim = c(nc, nsCh, p ^ 2))
    if (returnSe) {
      QInv <- array(dim = c(nsCh, (nc + 1) * p, (nc + 1) * p))
    }
    ## 'Fill' the SNP-dependent quantities, looping over the individuals.
    VBetaSnp <- X2[ch, ] %*% vInvArrRed
    for (i in 1:n) {
      VSnp <- VSnp + tcrossprod(X[ch, i], vInvArr[i, , ] %*% Y[, i])
    }
    for (cv in 1:nc) {
      X2VinvX1Arr[cv, , ] <- X[ch, ] %*% (vInvArrRed * as.numeric(W[cv, ]))
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
        QSnpInv <- solve(rbind(cbind(VBeta, X2VinvX1),
                               cbind(t(X2VinvX1), solve(VBetaSnpInv))))
        EffSe[, snpPos0 + snp] <- sqrt(diag(QSnpInv)[-(1:(p * nc))])
        QInv[snp, , ] <- QSnpInv
      }
    }
    if (returnSe) {
      VSnpQ <- matrix(data = 0, nrow = nsCh, ncol = p)
      for (i in 1:n) {
        VSnpQ <- VSnpQ + X[ch, i] %o% YSnpVInv[, i]
      }
      for (snp in 1:nsCh) {
        QSnp <- c(VQ, VSnpQ[snp, ])
        SS1[snpPos0 + snp] <- qScal - crossprod(QSnp, QInv[snp, , ] %*% QSnp)
      }
    }
    if (estCom) {
      ## Define SNP-dependent quantities
      X2VinvX1ArrCom <- array(data = 0, dim = c(nc, nsCh, p))
      if (returnSe) {
        QInvCom <- array(dim = c(nsCh, nc * p + 1, nc * p + 1))
      }
      ## 'Fill' the SNP-dependent quantities, looping over the individuals.
      VSnpCom <- X[ch, ] %*% s1
      VBetaSnpCom <- X2[ch, ] %*% s2
      for (cv in 1:nc) {
        X2VinvX1ArrCom[cv, , ] <- X[ch, ] %*% (t(S3) * as.numeric(W[cv, ]))
      }
      for (snp in 1:nsCh) {
        VBetaSnpInvCom <- 1 / VBetaSnpCom[snp]
        X2VinvX1Com <- matrix(data = t(X2VinvX1ArrCom[, snp, ]), ncol = nc * p,
                              byrow = TRUE)
        XSInvCom <- solve(VBeta - VBetaSnpInvCom * crossprod(X2VinvX1Com))
        EffCovComSnp <- XSInvCom %*%
          (v - VBetaSnpInvCom * VSnpCom[snp] * t(X2VinvX1Com))
        EffCom[snpPos0 + snp] <- VBetaSnpInvCom %*%
          (VSnpCom[snp] - X2VinvX1Com %*% EffCovComSnp)
        if (returnSe) {
          QSnpInvCom <- solve(rbind(cbind(VBeta, t(X2VinvX1Com)),
                                    c(X2VinvX1Com, 1 / VBetaSnpInvCom)))
          EffSeCom[snpPos0 + snp] <- sqrt(diag(QSnpInvCom)[-(1:(p * nc))])
          QInvCom[snp, , ] <- QSnpInvCom
        }
      }
      if (returnSe) {
        vSnpQCom <- rowSums(X[ch, ] %*% colSums(YSnpVInv))
        for (snp in 1:nsCh) {
          QSnpCom <- c(VQ, vSnpQCom[snp])
          SS1Com[snpPos0 + snp] <- qScal -
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
      FValsCom <- (SS0 - SS1Com) / SS1Com * dfCom
      pValsCom <- pf(q = FValsCom, df1 = 1, df2 = dfCom, lower.tail = FALSE)
      FValsQtlE <- ((SS1Com - SS1) / SS1) * dfFull / (p - 1)
      pValsQtlE <- pf(q = FValsQtlE, df1 = p - 1, df2 = dfFull,
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
      out$pValsCom <- pValsCom
      out$PValsQtlE <- pValsQtlE
    }
  }
  return(out)
}
