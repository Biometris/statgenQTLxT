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
#' param return.SE (Default TRUE) Should standard errors and p-values be returned ?
#' param  estimate.common (Default FALSE) Determines if also the common SNP-effect model is fitted.
#                          Under development: p-values and standard errors not yet implemented.

estimateEffects_nw <- function(Y,
                               W = matrix(data = 1, nrow = nrow(Y)),
                               X,
                               Vg,
                               Ve,
                               K,
                               return.SE = TRUE,
                               estimate.common = FALSE) {
  Y <- as.matrix(Y)
  W <- as.matrix(W)
  X <- as.matrix(X)
  # number of traits
  p <- ncol(Y)
  # number of individuals
  n <- nrow(Y)
  # number of covariates, including the intercept
  nc <- ncol(W)
  # number of SNPs
  ns <- ncol(X)
  ## Extract names of SNPs and individuals.
  snpNames <- colnames(X)
  genoNames <- colnames(Y)
  ############
  # eigen-decomposition of K
  w <- eigen(K, symmetric = TRUE)
  Dk <- diag(w$values)
  Uk <- w$vectors
  #################
  # transform Y, W and X
  Y <- t(Y) %*% Uk
  # transform covariates
  W <- t(W) %*% Uk
  # transform marker matrix
  X <- t(X) %*% Uk
  # square each element of X
  X2 <- X ^ 2
  ####################
  # Define V.inv.array
  V.inv.array <- make.V.inv.array(Vg = Vg, Ve = Ve , Dk = Dk)
  # create a second instance of V.inv.array, in matrix form
  V.inv.array.red <- V.inv.array
  dim(V.inv.array.red) <- c(n, p ^ 2)
  ###########################
  # Compute quantities that are independent of the SNPs
  Vbeta <- matrix(rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(tcrossprod(W[, m]), V.inv.array[m, , ])
  })), ncol = p * nc)
  v <- rowSums(sapply(X = 1:n, FUN = function(m) {
    kronecker(W[, m], V.inv.array[m, , ] %*% Y[, m])
  }))
  #####################################################
  # Define SNP-dependent quantities
  v.snp <- matrix(data = 0, nrow = ns, ncol = p)
  Vbeta.snp <- matrix(data = 0, nrow = ns, ncol = p ^ 2)
  X2VinvX1.array <- array(data = 0, dim = c(nc, ns, p ^ 2))
  EFFECTS <- matrix(data = 0, nrow = p, ncol = ns,
                    dimnames = list(genoNames, snpNames))
  if (return.SE) {
    EFFECTS.se <- EFFECTS
    EFFECTS.cov <- matrix(data = 0, nrow = nc * p, ncol = ns)
    Q.matrix.inv.all <- matrix(data = 0, nrow = ((nc + 1) * p) ^ 2, ncol = ns)
  }
  ################################################
  # 'Fill' the SNP-dependent quantities, looping over the individuals
  for (i in 1:n) {
    v.snp <- v.snp + tcrossprod(X[, i], V.inv.array[i, , ] %*% Y[, i])
    Vbeta.snp <- Vbeta.snp + tcrossprod(X2[, i], V.inv.array.red[i, ])
    for (cv in 1:nc) {
      X2VinvX1.array[cv, , ] <- X2VinvX1.array[cv, , ] +
        W[cv, i] * tcrossprod(X[, i], V.inv.array.red[i, ])
    }
  }
  ############################################
  # Now do the remaining calculations, looping over the SNPs
  for (snp in 1:ns) {
    Vbeta.snp.inv <- solve(matrix(data = Vbeta.snp[snp, ], ncol = p))
    X2VinvX1 <- matrix(t(X2VinvX1.array[, snp, ]), ncol = p, byrow = TRUE)
    X.S.inv <- solve(Vbeta - X2VinvX1 %*% Vbeta.snp.inv %*% t(X2VinvX1))
    EFFECTS[, snp] <- Vbeta.snp.inv %*%
      (v.snp[snp, ] + t(X2VinvX1) %*% X.S.inv %*%
         (X2VinvX1 %*% Vbeta.snp.inv %*% v.snp[snp, ] - v))
    if (return.SE) {
      EFFECTS.cov[, snp] <- X.S.inv %*% (v - X2VinvX1 %*% Vbeta.snp.inv %*%
                                           v.snp[snp, ])
      Q.matrix.inv <- solve(rbind(cbind(Vbeta, X2VinvX1),
                                  cbind(t(X2VinvX1), solve(Vbeta.snp.inv))))
      EFFECTS.se[, snp] <- sqrt(diag(Q.matrix.inv)[-(1:(p * nc))])
      Q.matrix.inv.all[, snp] <- as.numeric(Q.matrix.inv)
    }
  }
  ##############################
  if (return.SE) {
    pvalues <- setNames(rep(x = 1, times = ns), snpNames)
    df.full <- (n - (nc + 1)) * p
    df.red  <- (n - nc) * p
    ### Null model with the trait specific means only
    # (which should be the model for which the Vg and Ve estimates were obtained)
    est0 <- solve(Vbeta, v)
    fitted.mean0 <- matrix(est0, ncol = length(est0) / p) %*% W
    SS0 <- LL.quad.form.diag(Y = Y - fitted.mean0, V.inv.array = V.inv.array)
    q.scal <- rep(x = 0, times = ns)
    v.q <- matrix(data = 0, nrow = ns, ncol = p * nc)
    v.snp.q <- matrix(data = 0, nrow = ns, ncol = p)
    for (i in 1:n) {
      res.i <- t(matrix(data = 1, nrow = p, ncol = ns) * Y[, i]) -
        X[, i] * t(EFFECTS) - W[1, i] * t(EFFECTS.cov[1:p, ])
      if (nc > 1) {
        for (extra.cov in 2:nc) {
          res.i <- res.i - W[extra.cov, i] * t(EFFECTS.cov[(extra.cov - 1) * p + 1:p, ])
        }
      }
      res.iV.inv.i <- res.i %*% V.inv.array[i, , ]
      q.scal <- q.scal + rowSums(res.i * res.iV.inv.i)
      v.snp.q <- v.snp.q + X[, i] * res.iV.inv.i
      v.q[, 1:p] <- v.q[, 1:p] + W[1, i] * res.iV.inv.i
      if (nc > 1) {
        for (extra.cov in 2:nc) {
          v.q[, (extra.cov - 1) * p + 1:p] <- v.q[, (extra.cov - 1) * p + 1:p] +
            W[extra.cov, i] * res.iV.inv.i
        }
      }
    }
    for (snp in 1:ns) {
      Q.matrix.inv <- matrix(Q.matrix.inv.all[, ns], ncol = p * (nc + 1))
      q.snp <- matrix(c(v.q[snp, ], v.snp.q[snp, ]))
      SS1 <- q.scal[snp] - t(q.snp) %*% Q.matrix.inv %*% q.snp
      Fstat <- ((SS0 - SS1) / SS1) * df.full / (df.red - df.full)
      pvalues[snp] <- pf(q = Fstat, df1 = df.red - df.full,
                         df2 = df.full, lower.tail = FALSE)
    }
  }
  ##########################
  if (estimate.common) {
    # Define SNP-dependent quantities
    v.snp.common <- Vbeta.snp.common <- matrix(data = 0, nrow = ns)
    X2VinvX1.array.common <- array(data = 0, dim = c(nc, ns, p))
    EFFECTS.common <- matrix(data = 0, ncol = ns)
    if (return.SE) {
      EFFECTS.common.se <- EFFECTS.common
      EFFECTS.cov <- matrix(data = 0, nrow = nc * p, ncol = ns)
      Q.matrix.inv.all.common <- matrix(data = 0, nrow = (nc * p + 1) ^ 2,
                                        ncol = ns)
    }
    ################################################
    # 'Fill' the SNP-dependent quantities, looping over the individuals
    s1 <- s2 <- numeric(n)
    s3 <- matrix(nrow = n, ncol = p)
    for (i in 1:n) {
      s1[i] <- sum(V.inv.array[i, , ] %*% Y[, i])
      s2[i] <- sum(V.inv.array[i, , ])
      s3[i, ] <- rowSums(V.inv.array[i, , ])
    }
    for (i in 1:n) {
      v.snp.common <- v.snp.common + s1[i] * X[, i]
      Vbeta.snp.common <- Vbeta.snp.common + s2[i] * X2[, i]
      for (cv in 1:nc) {
        X2VinvX1.array.common[cv, , ] <- X2VinvX1.array.common[cv, , ] +
          W[cv, i] * tcrossprod(X[, i], s3[i, ])
      }
    }
    for (snp in 1:ns) {
      Vbeta.snp.inv.common <- matrix(data = 1 / Vbeta.snp.common[snp, ])
      X2VinvX1.common <- matrix(data = X2VinvX1.array.common[, snp, ],
                                ncol = p, byrow = TRUE)
      X.S.inv.common <- solve(Vbeta - t(X2VinvX1.common) %*%
                                Vbeta.snp.inv.common %*% X2VinvX1.common)
      EFFECTS.common[, snp] <- Vbeta.snp.inv.common %*%
        (v.snp.common[snp, ] + X2VinvX1.common %*% X.S.inv.common %*%
           (t(X2VinvX1.common) %*% Vbeta.snp.inv.common %*% v.snp.common[snp, ] - v))
      if (return.SE) {
        EFFECTS.cov[, snp] <- X.S.inv.common %*%
          (v - t(X2VinvX1.common) %*% Vbeta.snp.inv.common %*%
             v.snp.common[snp, ])
        Q.matrix.inv.common <- solve(rbind(cbind(Vbeta, t(X2VinvX1.common)),
                                           cbind((X2VinvX1.common),
                                                 solve(Vbeta.snp.inv.common))))
        EFFECTS.common.se[, snp] <- sqrt(diag(Q.matrix.inv.common)[-(1:(p * nc))])
        Q.matrix.inv.all.common[, snp] <- as.numeric(Q.matrix.inv.common)
      }
    }
  }    # end if (estimate.common) {
  #########################
  # return results
  # EFFECTS.common may also be returned
  if (return.SE) {
    return(list(effect.estimates = EFFECTS, effect.se = EFFECTS.se, pvalues = pvalues))
  } else {
    return(list(effect.estimates = EFFECTS))
  }
}
