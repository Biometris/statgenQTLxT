# FA = factor analytic ; a variation on the more general function EM_function_LS (low rank + sparse)
#
EMFA <- function(Y,
  K,
  X = matrix(rep(1, nrow(K))),
  CmHet = FALSE,
  DmHet = FALSE,
  tolerance = 1e-4,
  maxIter = 300L,
  CmStart = NULL,
  DmStart = NULL,
  mG = 1,
  mE = 1,
  maxDiag = 1e4,
  prediction = TRUE,
  stopIfDecreasing = FALSE,
  computeLogLik = FALSE) {


  # CmStart=NULL; DmStart=NULL; mG=1; mE=1; CmHet=T; DmHet=T;computeLogLik=TRUE; stopIfDecreasing=FALSE

  # Y=Y;K=K;X=X;maxIter=maxIter;tolerance=tolerance;CmStart=NULL;DmStart=NULL;mG=mG;mE=mE;CmHet=TRUE;DmHet=TRUE;maxDiag=100; prediction=TRUE;stopIfDecreasing=FALSE;computeLogLik=T

  # X default used to be : X=data.frame()

  # Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
  #                without missing values; NOT transformed
  # K            : the n x n kinship matrix; NOT transformed
  # X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
  #                NOT transformed
  # Dm.is.diagonal:
  # tolerance       : tolerance in the EM-algorithm : stop when the increase in log-lik. is less than tol
  # maxIter  : maximum number of iterations in the EM algorithm
  # CmStart     : starting values for Cm, as p x p matrix
  # DmStart     : starting values for Dm, as p x p matrix

  # also check : missing data

  Y <- as.matrix(Y)

  stopifnot(nrow(Y) == nrow(K))
  stopifnot(ncol(K) == nrow(K))

  nc <- ncol(X)
  if (nc > 0) {stopifnot(nrow(X) == nrow(K))}

  n <- ncol(K)
  p <- ncol(Y)

  ## check if mG and mE have sensible values
  if (mG != round(mG)) {stop("mG needs to be integer")}
  if (mE != round(mE)) {stop("mE needs to be integer")}
  if (mG < 0) {stop("mG cannot be negative")}
  if (mE < 0) {stop("mE cannot be negative")}
  if (mG >= p) {stop("mG needs to be smaller than the number of traits or environments")}
  if (mE >= p) {stop("mE needs to be smaller than the number of traits or environments")}

  if (nc > 0) {
    B <- matrix(0, nc, p) # the c x p matrix of coefficients (p traits)
    XtXinvXt <- solve(t(X) %*% X) %*% t(X)
  } else {
    B <- NULL
  }

  R <- Ginv(K)
  w <- eigen(R, symmetric = TRUE)
  Uk <- w$vectors
  Dk <- diag(w$values)
  lambdaR <- as.matrix(Dk)

  ## Set start values
  if (is.null(CmStart)) {
    if (mG == 0) {
      Cm <- diag(x = 2, nrow = p)
    } else {
      Cm <- Ginv((cor(Y) + diag(p)) / 4)
    }
  } else {
    Cm <- CmStart
  }

  if (is.null(DmStart)) {
    if (mE == 0) {
      Dm <- diag(x = 2, nrow = p)
    } else {
      Dm <- Ginv((cor(Y) + diag(p)) / 4)
    }
  } else {
    Dm <- DmStart
  }

  ## the model is Cm^{-1} = P^{-1} + W W^t =
  ## Given a starting value for Cm, set starting values for P and W
  if (mG > 0) {
    eigenC <- eigen(Ginv(Cm), symmetric = TRUE)
    UG <- as.matrix(eigenC$vectors[, 1:mG])
    psiG <- mean(eigenC$values[-(1:mG)])
    if (mG > 1) {
      rootLambdaG <- matrixRoot(diag(eigenC$values[1:mG] - psiG))
    } else {
      rootLambdaG <- matrix(sqrt(eigenC$values[1:mG] - psiG))
    }
    WG <- UG %*% rootLambdaG
    PG <- diag(p) / psiG
  } else {
    WG <- NULL
    PG <- NULL
  }
  if (mE > 0) {
    eigenD <- eigen(Ginv(Dm), symmetric = TRUE)
    UE <- as.matrix(eigenD$vectors[, 1:mE])
    psiE <- mean(eigenD$values[-(1:mE)])
    if (mE > 1) {
      rootLambdaE <- matrixRoot(diag(eigenD$values[1:mE] - psiE))
    } else {
      rootLambdaE <- matrix(sqrt(eigenD$values[1:mE] - psiE))
    }
    WE <- UE %*% rootLambdaE
    PE <- diag(p) / psiE
  } else {
    WE <- NULL
    PE <- NULL
  }

  ############
  continue <- TRUE
  decreased <- FALSE
  iter <- 1
  ELogLikCm <- ELogLikDm <- ELogLik <- -Inf
  mu <- matrix(rep(0, n * p), ncol = p)

  #############################################
  # EM
  while (continue & iter < maxIter) {
    ## Prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- as.matrix(Cm + t(Cm)) / 2
    Dm <- as.matrix(Dm + t(Dm)) / 2

    DmSqrtInv <- matrixRoot(Ginv(Dm))
    w1 <- eigen(DmSqrtInv %*% Cm %*% DmSqrtInv, symmetric = TRUE)
    Q1 <- w1$vectors
    lambda1 <- w1$values

    CmSqrtInv <- matrixRoot(Ginv(Cm))
    w2 <- eigen(CmSqrtInv %*% Dm %*% CmSqrtInv, symmetric = TRUE)
    Q2 <- w2$vectors
    lambda2 <- w2$values

    ## In the preprint of Dahl et al (arxiv, version 6 dec. 2013),
    # part1-part4 are the quantities, on the bottom part of p.6, in this order
    # Each time we compute the right hand side of the equation
    # In dahl_etal_2013_debug.r we checked the left hand side(s) as well
    # Also S1 and S2 correspond to p. 6 of their preprint

    if (nc > 0) {
      tUYminXb <- t(Uk) %*% (Y - X %*% B)
      S1 <- vecInvDiag(x = lambda1, y = diag(lambdaR)) * (tUYminXb %*% (matrixRoot(Dm) %*% Q1))
      S2 <- vecInvDiag(x = lambda2, y = 1 / diag(lambdaR)) * (tUYminXb %*% (MatrixRoot(Cm) %*% Q2))
    } else {
      S1 <- vecInvDiag(x = lambda1, y = diag(lambdaR)) * (t(Uk) %*% Y %*% MatrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / diag(lambdaR)) * (t(Uk) %*% Y %*% MatrixRoot(Cm) %*% Q2)
    }

    trP1 <- tracePInvDiag(x = lambda1, y = diag(lambdaR))
    trP2 <- tracePInvDiag(x = lambda2, y = 1 / diag(lambdaR))
    if (p > 1) {
      part1 <- DmSqrtInv %*% Q1 %*% (diag(trP1) %*% t(Q1) %*% DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% (diag(trP2) %*% t(Q2) %*% CmSqrtInv)
    } else {
      part1 <- DmSqrtInv %*% Q1 %*% (matrix(trP1) %*% t(Q1) %*% DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% (matrix(trP2) %*% t(Q2) %*% CmSqrtInv)
    }
    part3 <- CmSqrtInv %*% Q2 %*% t(S2) %*% t(CmSqrtInv %*% (Q2 %*% t(S2)))
    part4 <- DmSqrtInv %*% Q1 %*% t(S1) %*% lambdaR %*% t(DmSqrtInv %*% (Q1 %*% t(S1)))

    if (nc > 0) {
      mu <- matrix(Uk %*% S1 %*% t(DmSqrtInv %*% Q1), ncol = p)
      B <- XtXinvXt %*% (Y - mu)
    }

    ################
    # Compare with the 'naive' expressions:
    #Sigma <- ginv(kronecker(Dm,diag(n)) + kronecker(Cm,R)); M <- matrix(Sigma %*% (kronecker(Dm,diag(n))) %*% matrix(as.numeric(Y)),ncol=p)
    #part1.check <- trace.p(Sigma,p,n); part2.check <-  trace.p(kronecker(diag(p),R) %*% Sigma,p,n)
    #part3.check <- t(Y - M) %*% (Y - M); part4.check <- t(M) %*% R %*% M
    #part1.check-part1;part2.check-part2;part3.check-part3;part4.check-part4
    #part1;part2;part3;part4
    ##############################

    Omega1 <- as.matrix((part1 + part3) / n)
    Omega2 <- as.matrix((part2 + part4) / n)

    ################  update C

    if (mG == 0) {
      # Recall that the model is Cm^{-1} = P^{-1} + W W^t.
      #
      # when mG==0, W=0 and Cm=P
      if (p > 1) {
        if (CmHet) {
          PgNew <- diag(pmin(maxDiag, 1 / diag(Omega2)))
        } else {
          tau <- min(maxDiag, p / sum(diag(Omega2)))
          PgNew <- tau * diag(p)
        }
      } else {
        PgNew <- matrix(1 / as.numeric(Omega2))
      }
      WgNew <- NULL
      CmNew  <- PgNew
    } else {
      # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- matrixRoot(Omega2)
      A <- A * sqrt(nrow(A))

      if (!CmHet) {
        CmNewOutput <- update_FA_homogeneous_var(S = Omega2, m = mG)
        CmNewOutput$P <- diag(p) / CmNewOutput$sigma2
      } else {
        CmNewOutput <- updateFA(Y=A,
          WStart= WG,
          PStart= PG,
          hetVar= CmHet,
          maxDiag= maxDiag)
      }
      WgNew <- CmNewOutput$W
      PgNew <- CmNewOutput$P
      CmNew <- Ginv(Ginv(PgNew) + WgNew %*% t(WgNew))
    }

    ################  update D

    if (mE == 0) {
      # Recall that the model is Dm^{-1} = P^{-1} + W W^t.
      #
      # when mE==0, W=0 and Cm=P
      if (p > 1) {
        if (DmHet) {
          PeNew <- diag(pmin(maxDiag, 1 / diag(Omega1)))
        } else {
          tau <- min(maxDiag, p / sum(diag(Omega1)))
          PeNew <- diag(x = tau, nrow = p)
        }
      } else {
        PeNew <- matrix(1 / as.numeric(Omega1))
      }
      WeNew <- NULL
      DmNew  <- PeNew
    } else {
      # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- matrixRoot(Omega1)
      A <- A * sqrt(nrow(A))
      if (!DmHet) {
        DmNew.output <- update_FA_homogeneous_var(S = Omega1, m = mE)
        DmNew.output$P <- diag(p) / DmNew.output$sigma2
      } else {
        DmNew.output <- updateFA(Y = A,
          WStart = WE,
          PStart = PE,
          hetVar = DmHet,
          maxDiag = maxDiag)
      }
      WeNew <- DmNew.output$W
      PeNew <- DmNew.output$P
      DmNew <- Ginv(Ginv(PeNew) + WeNew %*% t(WeNew))
    }

    #######################
    ELogLikOldCm <- ELogLikCm
    ELogLikOldDm <- ELogLikDm
    ELogLikOld <- ELogLikOldCm + ELogLikOldDm

    ELogLikCm <- n * log(det(Cm)) - n * sum(diag(Cm %*% Omega2))
    ELogLikDm <- n * log(det(Dm)) - n * sum(diag(Dm %*% Omega1))
    ELogLik <- ELogLikCm + ELogLikDm

    if (stopIfDecreasing & iter > 50) {
      if (ELogLik < ELogLikOld - 0.1) {
        continue <- FALSE
        decreased <- TRUE
      }
    }
    if (iter / 1000 == round(iter / 1000)) {
      CmDiff <- sum(abs(CmNew - Cm))
      DmDiff <- sum(abs(DmNew - Dm))
      cat('Iteration ', iter, ' : ', CmDiff, '  ', DmDiff, '    ', ELogLik,'\n')
    }

    #############

    Cm <- CmNew
    Dm <- DmNew
    WG <- WgNew
    WE <- WeNew
    PG <- PgNew
    PE <- PeNew
    iter <- iter + 1

    if (abs(ELogLik - ELogLikOld) < tolerance) {continue <- FALSE}
  }

  #################################

  if (computeLogLik) {
    VInvArray <- makeVInvArray(Vg = solve(Cm), Ve = solve(Dm), Dk = Dk)
    VArray <- makeVArray(Vg = solve(Cm), Ve = solve(Dm), Dk = Dk)
    if (nc > 0) {
      XTransformed <- t(X) %*% Uk
    } else {
      XTransformed <- data.frame()}
    logLik <- LL.diag(t(Y) %*% Uk, X = XTransformed, V.array = VArray, V.inv.array = VInvArray)
  } else {
    logLik <- NA
  }

  #################################

  # prevent that Cm, Dm become asymmetric because of numerical inaccuracies
  Cm <- (Cm + t(Cm)) / 2
  Dm <- (Dm + t(Dm)) / 2

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype', 1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait', 1:p)}

  if (prediction & nc == 0) {
    mu <- matrix(Uk %*% S1 %*% t(DmSqrtInv %*% Q1), ncol=p)
  }

  pred.frame <- data.frame(trait = rep(colnames(Y), each = n),
    genotype = rep(rownames(Y), p),
    predicted = as.numeric(mu))

  return(list(Cm = Cm, Dm = Dm, B = B, WG = WG, WE = WE,
    PG = PG, PE = PE, pred = pred.frame, logLik = ELogLik, log.lik2 = logLik,
    tolerance = tolerance, converged = (!continue),
    n = n, p = p, nc = nc, n.iter = iter, decreased = decreased))
}


