# To do : check data types (matrix versus data.frame)
#         Other checks, in particular if row/colnames match
#
# NOTE : X, X and Y : all UNtransformed
#
#' param Y : n x p matrix of phenotypic values, for n individuals and p traits
#' param W : n x c matrix of covariates, including an intercept. Default: only intercept
#' param X: ns x n matrix of marker scores. No missing data allowed; neither non-segregating markers
#     The columns should correspond to the rows of Y and W, and to K
#' param K : n x n genetic relatedness matrix
#' param Vg, Ve : p x p matrices with the genetic resp. environmental covariances
#          The row and column order should correspond to Y
#' param returnSe (Default TRUE) Should standard errors and p-values be returned ?
#' param  estCom (Default FALSE) Determines if also the common SNP-effect model is fitted.
#                          Under development: p-values and standard errors not yet implemented.

estimateEffects_nw <- function(Y,
                               W = matrix(data = 1, nrow = nrow(Y)),
                               X,
                               Vg,
                               Ve,
                               K,
                               returnSe = TRUE,
                               estCom = FALSE) {
  Y <- as.matrix(Y)
  W <- as.matrix(W)
  X <- as.matrix(X)
  ## Extract number of traits, individuals, covariates and SNPs
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
  vInvArr <- make.V.inv.array(Vg = Vg, Ve = Ve , Dk = Dk)
  ## create a second instance of vInvArr, in matrix form
  vInvArrRed <- vInvArr
  dim(vInvArrRed) <- c(n, p ^ 2)
  ## Compute quantities that are independent of the SNPs
  VBeta <- matrix(rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(tcrossprod(W[, m]), vInvArr[m, , ])
  })), ncol = p * nc)
  v <- rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(W[, m], vInvArr[m, , ] %*% Y[, m])
  }))
  ## Define SNP-dependent quantities
  VSnp <- matrix(data = 0, nrow = ns, ncol = p)
  VBetaSnp <- matrix(data = 0, nrow = ns, ncol = p ^ 2)
  X2VinvX1Arr <- array(data = 0, dim = c(nc, ns, p ^ 2))
  Eff <- matrix(data = 0, nrow = p, ncol = ns,
                dimnames = list(genoNames, snpNames))
  if (returnSe) {
    EffSe <- Eff
    EffCov <- matrix(data = 0, nrow = nc * p, ncol = ns)
    QInv <- matrix(data = 0, nrow = ((nc + 1) * p) ^ 2, ncol = ns)
  }
  ## 'Fill' the SNP-dependent quantities, looping over the individuals
  for (i in 1:n) {
    VSnp <- VSnp + tcrossprod(X[, i], vInvArr[i, , ] %*% Y[, i])
    VBetaSnp <- VBetaSnp + tcrossprod(X2[, i], vInvArrRed[i, ])
    for (cv in 1:nc) {
      X2VinvX1Arr[cv, , ] <- X2VinvX1Arr[cv, , ] +
        W[cv, i] * tcrossprod(X[, i], vInvArrRed[i, ])
    }
  }
  ## Now do the remaining calculations, looping over the SNPs
  for (snp in 1:ns) {
    VBetaSnpInv <- solve(matrix(data = VBetaSnp[snp, ], ncol = p))
    X2VinvX1 <- matrix(t(X2VinvX1Arr[, snp, ]), ncol = p, byrow = TRUE)
    XSInv <- solve(VBeta - X2VinvX1 %*% VBetaSnpInv %*% t(X2VinvX1))
    Eff[, snp] <- VBetaSnpInv %*%
      (VSnp[snp, ] + t(X2VinvX1) %*% XSInv %*%
         (X2VinvX1 %*% VBetaSnpInv %*% VSnp[snp, ] - v))
    if (returnSe) {
      EffCov[, snp] <- XSInv %*% (v - X2VinvX1 %*% VBetaSnpInv %*% VSnp[snp, ])
      QSnpInv <- solve(rbind(cbind(VBeta, X2VinvX1),
                             cbind(t(X2VinvX1), solve(VBetaSnpInv))))
      ## FOR DEBUGGING NAN WARNING
      if (any(diag(QSnpInv)[-(1:(p * nc))] < 0)) browser()
      ##
      EffSe[, snp] <- sqrt(diag(QSnpInv)[-(1:(p * nc))])
      QInv[, snp] <- as.numeric(QSnpInv)
    }
  }
  if (returnSe) {
    pVals <- setNames(rep(x = 1, times = ns), snpNames)
    dfFull <- (n - (nc + 1)) * p
    dfRed  <- (n - nc) * p
    ## Null model with the trait specific means only
    ## (which should be the model for which the Vg and Ve estimates were obtained)
    est0 <- solve(VBeta, v)
    fitMean0 <- matrix(est0, ncol = length(est0) / p) %*% W
    SS0 <- LL.quad.form.diag(Y = Y - fitMean0, V.inv.array = vInvArr)
    qScal <- rep(x = 0, times = ns)
    VQ <- matrix(data = 0, nrow = ns, ncol = p * nc)
    VSnpQ <- matrix(data = 0, nrow = ns, ncol = p)
    for (i in 1:n) {
      Res <- t(matrix(data = 1, nrow = p, ncol = ns) * Y[, i]) -
        X[, i] * t(Eff) - W[1, i] * t(EffCov[1:p, ])
      if (nc > 1) {
        for (ec in 2:nc) {
          Res <- Res - W[ec, i] * t(EffCov[(ec - 1) * p + 1:p, ])
        }
      }
      ResVInv <- Res %*% vInvArr[i, , ]
      qScal <- qScal + rowSums(Res * ResVInv)
      VSnpQ <- VSnpQ + X[, i] * ResVInv
      VQ[, 1:p] <- VQ[, 1:p] + W[1, i] * ResVInv
      if (nc > 1) {
        for (ec in 2:nc) {
          VQ[, (ec - 1) * p + 1:p] <- VQ[, (ec - 1) * p + 1:p] +
            W[ec, i] * ResVInv
        }
      }
    }
    for (snp in 1:ns) {
      QSnpInv <- matrix(QInv[, ns], ncol = p * (nc + 1))
      QSnp <- matrix(c(VQ[snp, ], VSnpQ[snp, ]))
      SS1 <- qScal[snp] - t(QSnp) %*% QSnpInv %*% QSnp
      FVal <- ((SS0 - SS1) / SS1) * dfFull / (dfRed - dfFull)
      pVals[snp] <- pf(q = FVal, df1 = dfRed - dfFull, df2 = dfFull,
                       lower.tail = FALSE)
    }
  }
  if (estCom) {
    ## Define SNP-dependent quantities
    VSnpCom <- VBetaSnpCom <- matrix(data = 0, nrow = ns)
    X2VinvX1ArrCom <- array(data = 0, dim = c(nc, ns, p))
    EffCom <- matrix(data = 0, ncol = ns)
    if (returnSe) {
      EffSeCom <- EffCom
      EffCovCom <- matrix(data = 0, nrow = nc * p, ncol = ns)
      QInvCom <- matrix(data = 0, nrow = (nc * p + 1) ^ 2, ncol = ns)
    }
    ## 'Fill' the SNP-dependent quantities, looping over the individuals
    s1 <- s2 <- numeric(n)
    S3 <- matrix(nrow = n, ncol = p)
    for (i in 1:n) {
      s1[i] <- sum(vInvArr[i, , ] %*% Y[, i])
      s2[i] <- sum(vInvArr[i, , ])
      S3[i, ] <- rowSums(vInvArr[i, , ])
    }
    for (i in 1:n) {
      VSnpCom <- VSnpCom + s1[i] * X[, i]
      VBetaSnpCom <- VBetaSnpCom + s2[i] * X2[, i]
      for (cv in 1:nc) {
        X2VinvX1ArrCom[cv, , ] <- X2VinvX1ArrCom[cv, , ] +
          W[cv, i] * tcrossprod(X[, i], S3[i, ])
      }
    }
    for (snp in 1:ns) {
      VBetaSnpInvCom <- matrix(data = 1 / VBetaSnpCom[snp, ])
      X2VinvX1Com <- matrix(data = X2VinvX1ArrCom[, snp, ], ncol = p,
                            byrow = TRUE)
      XSInvCom <- solve(VBeta - t(X2VinvX1Com) %*% VBetaSnpInvCom %*%
                          X2VinvX1Com)
      EffCom[, snp] <- VBetaSnpInvCom %*%
        (VSnpCom[snp, ] + X2VinvX1Com %*% XSInvCom %*%
           (t(X2VinvX1Com) %*% VBetaSnpInvCom %*% VSnpCom[snp, ] - v))
      if (returnSe) {
        EffCovCom[, snp] <- XSInvCom %*%
          (v - t(X2VinvX1Com) %*% VBetaSnpInvCom %*% VSnpCom[snp, ])
        QSnpInvCom <- solve(rbind(cbind(VBeta, t(X2VinvX1Com)),
                                  cbind((X2VinvX1Com), solve(VBetaSnpInvCom))))
        EffSeCom[, snp] <- sqrt(diag(QSnpInvCom)[-(1:(p * nc))])
        QInvCom[, snp] <- as.numeric(QSnpInvCom)
      }
    }
  }
  ## Construct output.
  out <- list(effects = Eff)
  if (returnSe) {
    out$effectsSe = EffSe
    out$pVals = pVals
  }
  if (estCom) {
    out$effectsCom = EffCom
    if (returnSe) {
      out$effectsSe = EffSeCom
    }
  }
  return(out)
}
