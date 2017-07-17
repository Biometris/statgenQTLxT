#' Factor analytic variation of EM algoritm
#'
#' Implementation of the factor analytic variation of the EM algoritm as proposed by Dahl et al.
#'
#' @param Y an n x p matrix of observed phenotypes, on p traits or environments for n individuals.
#' No missing values are allowed.
#' @param K an n x n kinship matrix.
#' @param X an n x c covariate matrix, c being the number of covariates and n being the number
#' of genotypes. c has to be at least one (typically an intercept). No missing values are allowed.
#' If not provided a vector of 1s is used.
#' @param CmHet should an extra diagonal part is added in the model for the
#' precision matrix Cm?
#' @param DmHet should an extra diagonal part is added in the model for the
#' precision matrix Dm?
#' @param tolerance a numerical value. The iterating process stops if the difference in conditional
#' log-likelihood between two consecutive iterations drops below tolerance.
#' @param maxIter a numerical value for the maximum number of iterations.
#' @param CmStart a p x p matrix containing starting values for the precision matrix Cm.
#' @param DmStart a p x p matrix containing starting values for the precision matrix Dm.
#' @param mG an integer. The order of the genetic part of the model.
#' @param mE an integer. The order of the environmental part of the model.
#' @param maxDiag a numical value. The maximal value of the diagonal elements in the precision matrices
#' Cm and Dm (ignoring the low-rank part W W^t)
#' @param prediction should predicted values for Y be returned?
#' @param stopIfDecreasing should the iterating process stop if after 50 iterations the
#' log-likelihood decreases between two consecutive iterations?
#' @param computeLogLik should the log-likelihood be returned?
#'
#' @return A list containing results of the algoritm.
#'
#' @references Dahl et al. (2014). Network inference in matrix-variate Gaussian models with
#' non-indenpent noise.

## TO DO: also check : missing data
## TO DO: describe output. Everything needed?

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

  if (!is.matrix(Y)) Y <- as.matrix(Y)
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
    XtXinvXt <- solve(t(X) %*% X, t(X))
  } else {
    B <- NULL
  }

  w <- eigen(solve(as(K, "symmetricMatrix")), symmetric = TRUE)
  Uk <- w$vectors
  Dk <- diag(w$values)
  lambdaR <- as.matrix(Dk)

  ## Set starting values for Cm
  if (is.null(CmStart)) {
    if (mG == 0) {
      Cm <- as(diag(x = 2, nrow = p), "symmetricMatrix")
    } else {
      Cm <- solve(as((cor(Y) + diag(p)) / 4, "symmetricMatrix"))
    }
  } else {
    Cm <- CmStart
  }

  ## Set starting values for Dm
  if (is.null(DmStart)) {
    if (mE == 0) {
      Dm <- as(diag(x = 2, nrow = p), "symmetricMatrix")
    } else {
      Dm <- solve(as((cor(Y) + diag(p)) / 4, "symmetricMatrix"))
    }
  } else {
    Dm <- DmStart
  }
  ## The model is Cm^{-1} = P^{-1} + W W^t
  ## Given a starting value for Cm, set starting values for P and W
  if (mG > 0) {
    eigenC <- eigen(solve(Cm), symmetric = TRUE)
    Ug <- as.matrix(eigenC$vectors[, 1:mG])
    psiG <- mean(eigenC$values[-(1:mG)])
    if (mG > 1) {
      rootLambdaG <- matrixRoot(diag(eigenC$values[1:mG] - psiG))
    } else {
      rootLambdaG <- matrix(sqrt(eigenC$values[1:mG] - psiG))
    }
    Wg <- Ug %*% rootLambdaG
    Pg <- diag(x = 1 / psiG, nrow = p)
  } else {
    Wg <- NULL
    Pg <- NULL
  }
  ## The model is Dm^{-1} = P^{-1} + W W^t
  ## Given a starting value for Dm, set starting values for P and W
  if (mE > 0) {
    eigenD <- eigen(solve(Dm), symmetric = TRUE)
    Ue <- as.matrix(eigenD$vectors[, 1:mE])
    psiE <- mean(eigenD$values[-(1:mE)])
    if (mE > 1) {
      rootLambdaE <- matrixRoot(diag(eigenD$values[1:mE] - psiE))
    } else {
      rootLambdaE <- matrix(sqrt(eigenD$values[1:mE] - psiE))
    }
    We <- Ue %*% rootLambdaE
    Pe <- diag(x = 1 / psiE, nrow = p)
  } else {
    We <- NULL
    Pe <- NULL
  }

  continue <- TRUE
  decreased <- FALSE
  iter <- 1
  ELogLikCm <- ELogLikDm <- -Inf
  mu <- matrix(rep(0, n * p), ncol = p)

  ## EM following Dahl et al.
  while (continue && iter < maxIter) {
    ## Prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- (Cm + t(Cm)) / 2
    Dm <- (Dm + t(Dm)) / 2

    DmSqrtInv <- matrixRoot(solve(Dm))
    w1 <- eigen(DmSqrtInv %*% Cm %*% DmSqrtInv, symmetric = TRUE)
    Q1 <- w1$vectors
    lambda1 <- w1$values

    CmSqrtInv <- matrixRoot(solve(Cm))
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
      S1 <- vecInvDiag(x = lambda1, y = w$values) * (tUYminXb %*% (matrixRoot(Dm) %*% Q1))
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) * (tUYminXb %*% (matrixRoot(Cm) %*% Q2))
    } else {
      S1 <- vecInvDiag(x = lambda1, y = w$values) * (t(Uk) %*% Y %*% matrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) * (t(Uk) %*% Y %*% matrixRoot(Cm) %*% Q2)
    }

    trP1 <- tracePInvDiag(x = lambda1, y = w$values)
    trP2 <- tracePInvDiag(x = lambda2, y = 1 / w$values)
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
    Omega1 <- as.matrix((part1 + part3) / n)
    Omega2 <- as.matrix((part2 + part4) / n)

    ################
    # Compare with the 'naive' expressions:
    #Sigma <- ginv(kronecker(Dm,diag(n)) + kronecker(Cm,R)); M <- matrix(Sigma %*% (kronecker(Dm,diag(n))) %*% matrix(as.numeric(Y)),ncol=p)
    #part1.check <- trace.p(Sigma,p,n); part2.check <-  trace.p(kronecker(diag(p),R) %*% Sigma,p,n)
    #part3.check <- t(Y - M) %*% (Y - M); part4.check <- t(M) %*% R %*% M
    #part1.check-part1;part2.check-part2;part3.check-part3;part4.check-part4
    #part1;part2;part3;part4
    ##############################

    ## Update C
    if (mG == 0) {
      ## Recall that the model is Cm^{-1} = P^{-1} + W W^t.
      ## when mG == 0, W = 0 and Cm = P
      if (p > 1) {
        if (CmHet) {
          PgNew <- diag(pmin(maxDiag, 1 / diag(Omega2)))
        } else {
          tau <- min(maxDiag, p / sum(diag(Omega2)))
          PgNew <- diag(x = tau, nrow = p)
        }
      } else {
        PgNew <- matrix(1 / as.numeric(Omega2))
      }
      WgNew <- NULL
      CmNew  <- PgNew
    } else {
      ## When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- matrixRoot(Omega2)
      A <- A * sqrt(nrow(A))
      if (!CmHet) {
        CmNewOutput <- updateFAHomVar(S = Omega2, m = mG)
      } else {
        CmNewOutput <- updateFA(Y = A,
          WStart= Wg,
          PStart= Pg,
          hetVar= CmHet,
          maxDiag= maxDiag)
      }
      WgNew <- CmNewOutput$W
      PgNew <- CmNewOutput$P
      CmNew <- Ginv(Ginv(PgNew) + WgNew %*% t(WgNew))
    }

    ## Update D
    if (mE == 0) {
      ## Recall that the model is Dm^{-1} = P^{-1} + W W^t.
      ## when mE == 0, W = 0 and Cm = P
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
      ## When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- matrixRoot(Omega1)
      A <- A * sqrt(nrow(A))
      if (!DmHet) {
        DmNewOutput <- updateFAHomVar(S = Omega1, m = mE)
      } else {
        DmNewOutput <- updateFA(Y = A,
          WStart = We,
          PStart = Pe,
          hetVar = DmHet,
          maxDiag = maxDiag)
      }
      WeNew <- DmNewOutput$W
      PeNew <- DmNewOutput$P
      DmNew <- Ginv(Ginv(PeNew) + WeNew %*% t(WeNew))
    }

    ## Compute log-likelihood and check stopping criteria
    ELogLikOld <- ELogLikCm + ELogLikDm
    ELogLikCm <- n * determinant(Cm)[[1]][1] - n * sum(diag(Cm %*% Omega2))
    ELogLikDm <- n * determinant(Dm)[[1]][1] - n * sum(diag(Dm %*% Omega1))
    ELogLik <- ELogLikCm + ELogLikDm

    if (stopIfDecreasing && iter > 50) {
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

    ## Update values for next iteration
    Cm <- CmNew
    Dm <- DmNew
    Wg <- WgNew
    We <- WeNew
    Pg <- PgNew
    Pe <- PeNew
    continue <- abs(ELogLik - ELogLikOld) >= tolerance && continue
    iter <- iter + 1
  }

  ## Compute log-likelihood
  if (computeLogLik) {
    VInvArray <- makeVInvArray(Vg <- solve(Cm), Ve <- solve(Dm), Dk = Dk)
    VArray <- makeVArray(Vg = Vg, Ve = Ve, Dk = Dk)
    if (nc > 0) {
      XTransformed <- t(X) %*% Uk
    } else {
      XTransformed <- data.frame()}
    logLik <- LLDiag(t(Y) %*% Uk, X = XTransformed, VArray = VArray, VInvArray = VInvArray)
  } else {
    logLik <- NA
  }

  ## prevent that Cm, Dm become asymmetric because of numerical inaccuracies
  Cm <- (Cm + t(Cm)) / 2
  Dm <- (Dm + t(Dm)) / 2

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype', 1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait', 1:p)}

  if (prediction & nc == 0) {
    mu <- matrix(Uk %*% S1 %*% t(DmSqrtInv %*% Q1), ncol = p)
  }

  predFrame <- data.frame(trait = rep(colnames(Y), each = n),
    genotype = rep(rownames(Y), p),
    predicted = as.numeric(mu))

  return(list(Cm = Cm, Dm = Dm, B = B, Wg = Wg, We = We,
    Pg = Pg, Pe = Pe, pred = predFrame, logLik = ELogLik, logLik2 = logLik,
    tolerance = tolerance, converged = (!continue),
    n = n, p = p, nc = nc, n.iter = iter, decreased = decreased))
}


