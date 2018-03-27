# To do : check data types (matrix versus data.frame)
#         Other checks, in particular if row/colnames match
#
# NOTE : X, mm and Y : all UNtransformed
#
# Y : n x p matrix of phenotypic values, for n individuals and p traits
# W : n x c matrix of covariates, including an intercept. Default: only intercept
# mm: ns x n matrix of marker scores. No missing data allowed; neither non-segregating markers
#     The columns should correspond to the rows of Y and W, and to K
# K : n x n genetic relatedness matrix
# Vg, Ve : p x p matrices with the genetic resp. environmental covariances
#          The row and column order should correspond to Y


estimateEffects_WK <- function(Y, W = matrix(1,nrow(Y),1), mm, Vg, Ve, K,
                               return.SE = TRUE, estimate.common = FALSE) {

  # number of traits
  p <- ncol(Y)

  # number of individuals
  n <- nrow(Y)

  # number of covariates, including the intercept
  nc   <- ncol(W)

  # number of SNPs
  ns <- nrow(mm)

  ############
  # eigen-decomposition of K

  w <- eigen(K)

  Dk <- diag(w$values)

  Uk <- w$vectors

  #################
  # transform Y, W and mm

  # transform phenotypes
  rownames.Y <- rownames(Y)
  Y <- t(Y) %*% Uk
  colnames(Y) <- rownames.Y

  # transform covariates
  W <- t(W) %*% Uk

  # transform marker matrix
  MM  <- as.matrix(mm) %*% Uk

  # square each element
  MM2 <- MM * MM

  ####################
  # Define V.inv.array

  V.inv.array <- make.V.inv.array(Vg=Vg, Ve=Ve , Dk=Dk)

  # create a second instance of V.inv.array, in matrix form
  V.inv.array.red <- V.inv.array
  dim(V.inv.array.red) <- c(n, p * p)

  ###########################
  # Compute quantities that are independent of the SNPs

  Vbeta <- matrix(apply(sapply(1:n,function(m){kronecker((matrix(W[,m])) %*% t(matrix(W[,m])),V.inv.array[m,,])}),1,sum),ncol=p*nc)

  v     <- matrix(apply(sapply(1:n,function(m){kronecker((matrix(W[,m])),V.inv.array[m,,] %*% matrix(Y[,m]))}),1,sum))

  #####################################################
  # Define SNP-dependent quantities

  v.snp <- matrix(0, ns, p)

  Vbeta.snp <- matrix(0, ns, p * p)

  X2VinvX1.array <- array(0, dim = c(nc, ns, p * p))

  EFFECTS <- matrix(0, p, ns)

  if (return.SE) {
    EFFECTS.cov <- matrix(0, nc * p, ns)
    EFFECTS.se <- matrix(0, p, ns)
    Q.matrix.inv.all <- matrix(0, ((nc+1) * p)^2, ns)
  }

  ################################################
  # 'Fill' the SNP-dependent quantities, looping over the individuals

  for (i in 1:n) {

    v.snp <- v.snp + tcrossprod(MM[,i], V.inv.array[i,,] %*% matrix(Y[,i]))

    Vbeta.snp <- Vbeta.snp + tcrossprod(MM2[,i], V.inv.array.red[i,])

    for (cv in 1:nc) {

      X2VinvX1.array[cv, , ] <- X2VinvX1.array[cv, , ] + W[cv, i] * tcrossprod(MM[,i], V.inv.array.red[i,])

    }

  }

  ############################################
  # Now do the remaining calculations, looping over the SNPs


  for (snp in 1:ns) {
      #snp=1
      Vbeta.snp.inv <- solve(matrix(Vbeta.snp[snp,], ncol=p))

      X2VinvX1 <- matrix(matrix(t(X2VinvX1.array[,snp,])), byrow=T, ncol=p)

      X.S.inv <- solve(Vbeta - X2VinvX1 %*% Vbeta.snp.inv %*% t(X2VinvX1))

      EFFECTS[, snp] <- Vbeta.snp.inv %*% matrix(v.snp[snp, ]) + Vbeta.snp.inv %*% t(X2VinvX1) %*% X.S.inv %*% (X2VinvX1 %*% Vbeta.snp.inv %*% matrix(v.snp[snp, ]) - v)

      if (return.SE) {
        EFFECTS.cov[, snp] <- X.S.inv %*% (v - X2VinvX1 %*% Vbeta.snp.inv %*% matrix(v.snp[snp, ]))
        Q.matrix.inv      <- solve(rbind(cbind(Vbeta, X2VinvX1), cbind(t(X2VinvX1), solve(Vbeta.snp.inv))))
        EFFECTS.se[, snp] <- sqrt(diag(Q.matrix.inv)[-(1:(p*(nc)))])
        Q.matrix.inv.all[, snp] <- as.numeric(Q.matrix.inv)
      }
  }


  ##############################

  if (return.SE) {

    pvalues <- rep(1, ns)

    df.full <- (n - (nc+1))*p

    df.red  <- (n - nc)*p

    ### Null model with the trait specific means only
    #    (which should be the model for which the Vg and Ve estimates were obtained)

    est0 <- as.numeric(solve(Vbeta) %*% v)

    fitted.mean0 <- matrix(est0, ncol= length(est0)/p ) %*% W

    SS0 <- LL.quad.form.diag(Y=Y-fitted.mean0, V.inv.array=V.inv.array)

    q.scal <- rep(0, ns)

    v.q <- matrix(0, ns, p * nc)

    v.snp.q <- matrix(0, ns, p)

    for (i in 1:n) {

      res.i <- t(matrix(1, p, ns) * Y[,i]) - (MM[,i] * t(EFFECTS)) - t(EFFECTS.cov[1:p,] * W[1,i])
      if (nc > 1) {
        for (extra.cov in 2:nc) {
          res.i <- res.i - t(EFFECTS.cov[(extra.cov-1)*p + 1:p,] * W[extra.cov,i])
        }
      }

      q.scal <- q.scal + apply(res.i * (res.i %*% V.inv.array[i,,]), 1, sum)
      v.snp.q <- v.snp.q + MM[,i] *(res.i %*% V.inv.array[i,,])

      v.q[,1:p] <- v.q[,1:p] + W[1,i] *(res.i %*% V.inv.array[i,,])
      if (nc > 1) {
        for (extra.cov in 2:nc) {
          v.q[,(extra.cov-1)*p + 1:p] <- v.q[,(extra.cov-1)*p + 1:p] + W[extra.cov,i] *(res.i %*% V.inv.array[i,,])
        }
      }
    }

    SS1.full <- rep(NA, ns)   # these values are stored, for later use for QTL x E

    for (snp in 1:ns) {

      Q.matrix.inv <- matrix(Q.matrix.inv.all[,ns], ncol = p * (nc + 1))

      q.snp <- matrix(c(v.q[snp,], v.snp.q[snp,]))

      SS1.full[snp] <- q.scal[snp] - t(q.snp) %*% Q.matrix.inv %*% q.snp

      Fstat <- ((SS0-SS1.full[snp]) / (SS1.full[snp])) * (df.full) / (df.red-df.full)

      pvalues[snp] <- pf(q=Fstat, df1=df.red-df.full, df2=df.full, lower.tail = F)
    }

  }

  ##########################

  if (estimate.common) {

    # Define SNP-dependent quantities

    v.snp.common <- matrix(0, ns, 1)

    Vbeta.snp.common <- matrix(0, ns, 1)

    X2VinvX1.array.common <- array(0, dim = c(nc, ns, p))

    EFFECTS.common <- matrix(0, 1, ns)

    if (return.SE) {
      EFFECTS.common.se <- matrix(0, 1, ns)
      EFFECTS.cov <- matrix(0, nc * p, ns)
      Q.matrix.inv.all.common <- matrix(0, (nc*p + 1)^2, ns)
    }

    ################################################
    # 'Fill' the SNP-dependent quantities, looping over the individuals

    s1 <- rep(0, n)
    s2 <- rep(0, n)
    s3 <- matrix(0, n, p)

    for (i in 1:n) {
      s1[i]  <- sum(V.inv.array[i,,] %*% matrix(Y[,i]))
      s2[i]  <- sum(V.inv.array[i,,])
      s3[i,] <- apply(V.inv.array[i,,],1,sum)
    }

    for (i in 1:n) {

      v.snp.common <- v.snp.common + s1[i] * MM[,i]

      Vbeta.snp.common <- Vbeta.snp.common + s2[i] * MM2[,i]

      for (cv in 1:nc) {
        X2VinvX1.array.common[cv, , ] <- X2VinvX1.array.common[cv, , ] + W[cv, i] * tcrossprod(MM[,i], s3[i,])
      }
    }

    for (snp in 1:ns) {

        Vbeta.snp.inv.common <- matrix(1/Vbeta.snp.common[snp,])

        X2VinvX1.common <- matrix(matrix(t(X2VinvX1.array.common[,snp,])), byrow=T, ncol=p)

        X.S.inv.common <- solve(Vbeta - t(X2VinvX1.common) %*% Vbeta.snp.inv.common %*% (X2VinvX1.common))

        EFFECTS.common[, snp] <- Vbeta.snp.inv.common %*% matrix(v.snp.common[snp, ]) + Vbeta.snp.inv.common %*% (X2VinvX1.common) %*% X.S.inv.common %*% (t(X2VinvX1.common) %*% Vbeta.snp.inv.common %*% matrix(v.snp.common[snp, ]) - v)

        if (return.SE) {
          EFFECTS.cov[, snp] <- X.S.inv.common %*% (v - t(X2VinvX1.common) %*% Vbeta.snp.inv.common %*% matrix(v.snp.common[snp, ]))
          Q.matrix.inv.common      <- solve(rbind(cbind(Vbeta, t(X2VinvX1.common)), cbind((X2VinvX1.common), solve(Vbeta.snp.inv.common))))
          EFFECTS.common.se[, snp] <- sqrt(diag(Q.matrix.inv.common)[-(1:(p*(nc)))])
          Q.matrix.inv.all.common[, snp] <- as.numeric(Q.matrix.inv.common)
        }
    }

    if (return.SE) {

      pvalues.common <- pvalues.qtlE <- rep(1, ns)

      df.full <- (n - (nc+1))*p

      df.red  <- (n - nc)*p

      df.common  <- df.red - 1

      ### Null model with the trait specific means only
      #    (which should be the model for which the Vg and Ve estimates were obtained)

      q.scal <- rep(0, ns)

      v.q <- matrix(0, ns, p * nc)

      v.snp.q <- matrix(0, ns, 1)

      for (i in 1:n) {

        res.i <- t(matrix(1, p, ns) * Y[,i]) - (matrix(MM[,i]) %*% matrix(1, 1, p)) * as.numeric(EFFECTS.common) - t(EFFECTS.cov[1:p,] * W[1,i])
        if (nc > 1) {
          for (extra.cov in 2:nc) {
            res.i <- res.i - t(EFFECTS.cov[(extra.cov-1)*p + 1:p,] * W[extra.cov,i])
          }
        }

        q.scal <- q.scal + apply(res.i * (res.i %*% V.inv.array[i,,]), 1, sum)

        v.snp.q <- v.snp.q + apply(res.i %*% V.inv.array[i,,], 1, sum) * as.numeric(MM[,i])

        v.q[,1:p] <- v.q[,1:p] + W[1,i] *(res.i %*% V.inv.array[i,,])
        if (nc > 1) {
          for (extra.cov in 2:nc) {
            v.q[,(extra.cov-1)*p + 1:p] <- v.q[,(extra.cov-1)*p + 1:p] + W[extra.cov,i] *(res.i %*% V.inv.array[i,,])
          }
        }
      }

      SS1.common <- rep(NA, ns)

      for (snp in 1:ns) {

        Q.matrix.inv <- matrix(Q.matrix.inv.all.common[,ns], ncol = p * nc + 1)

        q.snp <- matrix(c(v.q[snp,], v.snp.q[snp,]))

        SS1.common[snp] <- q.scal[snp] - t(q.snp) %*% Q.matrix.inv %*% q.snp

        Fstat <- ((SS0-SS1.common[snp]) / SS1.common[snp]) * (df.common) / (df.red-df.common)

        pvalues.common[snp] <- pf(q=Fstat, df1=df.red-df.common, df2=df.common, lower.tail = F)

        ###########

        Fstat.qtlE <- ((SS1.common[snp]-SS1.full[snp]) / (SS1.full[snp])) * (df.full) / (df.common-df.full)

        pvalues.qtlE[snp] <- pf(q=Fstat.qtlE, df1=df.common-df.full, df2=df.full, lower.tail = F)

      }

    }


  }    # end if (estimate.common) {

  #########################
  # return results

  if (estimate.common) {
    if (return.SE) {
      return(list(effect.estimates = EFFECTS, effect.se = EFFECTS.se,
                  pvalues = pvalues, effect.estimates.common = EFFECTS.common,
                  pvalues.common = pvalues.common, pvalues.qtlE = pvalues.qtlE
                  ))
                  # SS1.common = SS1.common, SS1.full = SS1.full, SS0 = SS0
    } else {
      return(list(effect.estimates = EFFECTS, effect.estimates.common = EFFECTS.common))
    }
  } else {
    if (return.SE) {
      return(list(effect.estimates = EFFECTS, effect.se = EFFECTS.se, pvalues = pvalues))
    } else {
      return(list(effect.estimates = EFFECTS))
    }
  }
}
