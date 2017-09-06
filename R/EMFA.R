#' Factor analytic variation of EM algoritm
#'
#' Implementation of the factor analytic variation of the EM algoritm as proposed by Dahl et al. (2013).
#'
#' @param Y an n x p matrix of observed phenotypes, on p traits or environments for n individuals.
#' No missing values are allowed.
#' @param K an n x n kinship matrix.
#' @param X an n x c covariate matrix, c being the number of covariates and n being the number
#' of genotypes. c has to be at least one (typically an intercept). No missing values are allowed.
#' If not provided a vector of 1s is used.
#' @param CmHet should an extra diagonal part be added in the model for the
#' precision matrix Cm?
#' @param DmHet should an extra diagonal part be added in the model for the
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
#' @return A list containing the following components
#' \itemize{
#' \item{\code{Cm} final value for the precision matrix Cm.}
#' \item{\code{Dm} final value for the precision matrix Dm.}
#' \item{\code{logLik} log-likelihood}
#' \item{\code{logLik2} log-likelihood as in Zhou and Stephens (2014)}
#' \item{\code{nIter} the number of iterations.}
#' \item{\code{converged} did the algorithm converge?}
#' \item{\code{decreased} did the algorithm stop because the log-likelihood decreased
#' between iterations.}
#' }
#'
#' @references Dahl et al. (2013). Network inference in matrix-variate Gaussian models with
#' non-independent noise. arXiv preprint arXiv:1312.1622.
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear mixed model algorithms for
#' genome-wide association studies. Nature Methods, February 2014, Vol. 11, p. 407–409
#'
#' @importFrom methods as
#'
#' @keywords internal

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
  ## Check input.
  if (missing(Y) || !is.matrix(Y) || anyNA(Y))
    stop("Y should be a matrix without missing values.")
  if (missing(K) || !is.matrix(K) || nrow(K) != nrow(Y) || ncol(K) != nrow(Y) || anyNA(K))
    stop("K should be a matrix without missing values with the same number of rows
      and columns as the number of rows in Y.")
  if(!is.matrix(X) || anyNA(X))
    stop("X should be a matrix without missing values.")
  nc <- ncol(X)
  if (nc > 0 && nrow(X) != nrow(K))
    stop("X and K should have the same number of rows.")
  n <- ncol(K)
  p <- ncol(Y)
  if (!is.numeric(mG) || mG != round(mG) || mG < 0 || mG > p)
    stop("mG should be a positive integer between 0 and the number of traits.")
  if (!is.numeric(mE) || mE != round(mE) || mE < 0 || mE > p)
    stop("mE should be a positive integer between 0 and the number of traits.")
  if (nc > 0) {
    B <- matrix(0, nc, p) # the c x p matrix of coefficients (p traits)
    XtXinvXt <- solve(crossprod(X), t(X))
  } else {
    B <- NULL
  }
  w <- eigen(solve(as(K, "symmetricMatrix")), symmetric = TRUE)
  Uk <- w$vectors
  Dk <- w$values
  lambdaR <- diag(w$values)
  ## Set starting values for Cm
  if (is.null(CmStart)) {
    if (mG == 0) {
      Cm <- Matrix::Diagonal(n = p, x = 2)
    } else {
      Cm <- solve(as((cor(Y) + diag(p)) / 4, "symmetricMatrix"))
    }
  } else {
    Cm <- CmStart
  }
  ## Set starting values for Dm
  if (is.null(DmStart)) {
    if (mE == 0) {
      Dm <- Matrix::Diagonal(n = p, x = 2)
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
  ## Set starting values.
  continue <- TRUE
  decreased <- FALSE
  iter <- 1
  ELogLikCm <- ELogLikDm <- -Inf
  mu <- matrix(rep(0, n * p), ncol = p)
  ## EM following the notation of Dahl et al.
  while (continue && iter < maxIter) {
    DmSqrtInv <- matrixRoot(solve(Dm))
    w1 <- eigen(DmSqrtInv %*% Cm %*% DmSqrtInv, symmetric = TRUE)
    Q1 <- w1$vectors
    lambda1 <- w1$values
    CmSqrtInv <- matrixRoot(solve(Cm))
    w2 <- eigen(CmSqrtInv %*% Dm %*% CmSqrtInv, symmetric = TRUE)
    Q2 <- w2$vectors
    lambda2 <- w2$values
    if (nc > 0) {
      tUYminXb <- crossprod(Uk, Y - X %*% B)
      S1 <- vecInvDiag(x = lambda1, y = w$values) * (tUYminXb %*% matrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) * (tUYminXb %*% matrixRoot(Cm) %*% Q2)
    } else {
      S1 <- vecInvDiag(x = lambda1, y = w$values) * crossprod(Uk, Y %*% matrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) * crossprod(Uk, Y %*% matrixRoot(Cm) %*% Q2)
    }
    trP1 <- tracePInvDiag(x = lambda1, y = w$values)
    trP2 <- tracePInvDiag(x = lambda2, y = 1 / w$values)
    if (p > 1) {
      part1 <- DmSqrtInv %*% Q1 %*% diag(trP1) %*% crossprod(Q1, DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% diag(trP2) %*% crossprod(Q2, CmSqrtInv)
    } else {
      part1 <- DmSqrtInv %*% Q1 %*% matrix(trP1) %*% crossprod(Q1, DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% matrix(trP2) %*% crossprod(Q2, CmSqrtInv)
    }
    part3 <- tcrossprod(CmSqrtInv %*% Q2 %*% t(S2))
    part4 <- tcrossprod(DmSqrtInv %*% Q1, S1) %*% lambdaR %*% t(tcrossprod(DmSqrtInv %*% Q1, S1))
    if (nc > 0) {
      mu <- matrix(tcrossprod(Uk %*% S1, DmSqrtInv %*% Q1), ncol = p)
      B <- XtXinvXt %*% (Y - mu)
    }
    Omega1 <- Matrix::forceSymmetric(as.matrix((part1 + part3) / n))
    Omega2 <- Matrix::forceSymmetric(as.matrix((part2 + part4) / n))
    ## Update Cm
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
          WStart = Wg,
          PStart = Pg,
          hetVar = CmHet,
          maxDiag = maxDiag)
      }
      WgNew <- CmNewOutput$W
      PgNew <- CmNewOutput$P
      CmNew <- MASS::ginv(MASS::ginv(PgNew) + tcrossprod(WgNew))
    }
    ## Update Dm
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
      DmNew <- MASS::ginv(MASS::ginv(PeNew) + tcrossprod(WeNew))
    }
    ## Compute log-likelihood and check stopping criteria
    ELogLikOld <- ELogLikCm + ELogLikDm
    ELogLikCm <- n * (Matrix::determinant(Cm)[[1]][1] - sum(Matrix::diag(Cm %*% Omega2)))
    ELogLikDm <- n * (Matrix::determinant(Dm)[[1]][1] - sum(Matrix::diag(Dm %*% Omega1)))
    ELogLik <- ELogLikCm + ELogLikDm
    if (stopIfDecreasing && iter > 50) {
      if (ELogLik < ELogLikOld - 0.1) {
        continue <- FALSE
        decreased <- TRUE
      }
    }
    if (iter %% 1000 == 0) {
      CmDiff <- sum(abs(CmNew - Cm))
      DmDiff <- sum(abs(DmNew - Dm))
      cat("Iteration ", iter, " : ", CmDiff, "  ", DmDiff, "    ", ELogLik,"\n")
    }
    ## Update values for next iteration
    ## Prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- Matrix::forceSymmetric(CmNew)
    Dm <- Matrix::forceSymmetric(DmNew)
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
      XTransformed <- crossprod(X, Uk)
    } else {
      XTransformed <- data.frame()
    }
    logLik <- LLDiag(crossprod(Y, Uk), X = XTransformed, VArray = VArray, VInvArray = VInvArray)
  } else {
    logLik <- NA
  }
  ## Add default names if needed.
  if (is.null(rownames(Y))) {rownames(Y) <- paste0("genotype", 1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0("trait", 1:p)}
  if (prediction & nc == 0) {
    mu <- matrix(tcrossprod(Uk %*% S1, DmSqrtInv %*% Q1), ncol = p)
  }
  predFrame <- data.frame(trait = rep(colnames(Y), each = n),
    genotype = rep(rownames(Y), p),
    predicted = as.numeric(mu))

  return(list(Cm = Cm, Dm = Dm, logLik = ELogLik, logLik2 = logLik, nIter = iter,
    converged = (!continue), decreased = decreased))
}


