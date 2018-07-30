#' Factor analytic variation of EM algoritm
#'
#' Implementation of the factor analytic variation of the EM algoritm as
#' proposed by Dahl et al. (2013).
#'
#' @param Y An n x p matrix of observed phenotypes, on p traits or environments
#' for n individuals. No missing values are allowed.
#' @param K An n x n kinship matrix.
#' @param X An n x c covariate matrix, c being the number of covariates and n
#' being the number of genotypes. c has to be at least one (typically
#' an intercept). No missing values are allowed. If not provided a vector of 1s
#' is used.
#' @param CmHet Should an extra diagonal part be added in the model for the
#' precision matrix Cm?
#' @param DmHet Should an extra diagonal part be added in the model for the
#' precision matrix Dm?
#' @param tolerance A numerical value. The iterating process stops if the
#' difference in conditional log-likelihood between two consecutive iterations
#' drops below tolerance.
#' @param maxIter A numerical value for the maximum number of iterations.
#' @param CmStart A p x p matrix containing starting values for the precision
#' matrix Cm.
#' @param DmStart A p x p matrix containing starting values for the precision
#' matrix Dm.
#' @param mG An integer. The order of the genetic part of the model.
#' @param mE An integer. The order of the environmental part of the model.
#' @param maxDiag A numical value. The maximal value of the diagonal elements
#' in the precision matrices Cm and Dm (ignoring the low-rank part W W^t)
#' @param prediction Should predicted values for Y be returned?
#' @param stopIfDecreasing Should the iterating process stop if after 50
#' iterations the log-likelihood decreases between two consecutive iterations?
#' @param computeLogLik Should the log-likelihood be returned?
#'
#' @return A list containing the following components
#' \itemize{
#' \item{\code{Cm} final value for the precision matrix Cm.}
#' \item{\code{Dm} final value for the precision matrix Dm.}
#' \item{\code{logLik} log-likelihood}
#' \item{\code{logLik2} log-likelihood as in Zhou and Stephens (2014)}
#' \item{\code{nIter} the number of iterations.}
#' \item{\code{converged} did the algorithm converge?}
#' \item{\code{decreased} did the algorithm stop because the log-likelihood
#' decreased between iterations.}
#' }
#'
#' @references Dahl et al. (2013). Network inference in matrix-variate Gaussian
#' models with non-independent noise. arXiv preprint arXiv:1312.1622.
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @importFrom methods as
#'
#' @keywords internal
EMFA <- function(Y,
                 K,
                 X = Matrix::Matrix(rep(1, nrow(K))),
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
  if (missing(Y) || !(is.matrix(Y) || inherits(Y, "Matrix")) || anyNA(Y)) {
    stop("Y should be a matrix without missing values.\n")
  }
  if (missing(K) || !(is.matrix(K) || inherits(K, "Matrix")) ||
      nrow(K) != nrow(Y) || ncol(K) != nrow(Y) || anyNA(K)) {
    stop(paste("K should be a matrix without missing values with the same",
               "number of rows and columns as the number of rows in Y.\n"))
  }
  if (!(is.matrix(X) || inherits(X, "Matrix")) || anyNA(X)) {
    stop("X should be a matrix without missing values.\n")
  }
  nc <- ncol(X)
  if (nc > 0 && nrow(X) != nrow(K)) {
    stop("X and K should have the same number of rows.")
  }
  n <- ncol(K)
  p <- ncol(Y)
  if (!is.numeric(mG) || mG != round(mG) || mG < 0 || mG > p) {
    stop("mG should be a positive integer between 0 and the number of traits.")
  }
  if (!is.numeric(mE) || mE != round(mE) || mE < 0 || mE > p) {
    stop("mE should be a positive integer between 0 and the number of traits.")
  }
  if (nc > 0) {
    B <- Matrix::Matrix(0, nc, p) # the c x p matrix of coefficients (p traits)
    XtXinvXt <- Matrix::solve(Matrix::crossprod(X), Matrix::t(X))
  } else {
    B <- NULL
  }
  w <- eigen(Matrix::solve(K), symmetric = TRUE)
  Uk <- as(w$vectors, "dgeMatrix")
  Dk <- w$values
  lambdaR <- Matrix::Diagonal(x = w$values)
  ## Set starting values for Cm
  if (is.null(CmStart)) {
    if (mG == 0) {
      Cm <- Matrix::Diagonal(n = p, x = 2)
    } else {
      Cm <- Matrix::solve(cor(as.matrix(Y)) + Matrix::Diagonal(n = p) / 4)
    }
  } else {
    Cm <- CmStart
  }
  ## Set starting values for Dm
  if (is.null(DmStart)) {
    if (mE == 0) {
      Dm <- Matrix::Diagonal(n = p, x = 2)
    } else {
      Dm <- Matrix::solve(cor(as.matrix(Y)) + Matrix::Diagonal(n = p) / 4)
    }
  } else {
    Dm <- DmStart
  }
  ## The model is Cm^{-1} = P^{-1} + W W^t
  ## Given a starting value for Cm, set starting values for P and W
  if (mG > 0) {
    eigenC <- eigen(Matrix::solve(Cm), symmetric = TRUE)
    Ug <- as(eigenC$vectors[, 1:mG], "dgeMatrix")
    psiG <- mean(eigenC$values[-(1:mG)])
    if (mG > 1) {
      rootLambdaG <- matrixRoot(
        Matrix::Diagonal(x = eigenC$values[1:mG] - psiG))
    } else {
      rootLambdaG <- as(sqrt(eigenC$values[1:mG] - psiG), "dgeMatrix")
    }
    Wg <- Ug %*% rootLambdaG
    Pg <- Matrix::Diagonal(n = p, x = 1 / psiG)
  } else {
    Wg <- NULL
    Pg <- NULL
  }
  ## The model is Dm^{-1} = P^{-1} + W W^t
  ## Given a starting value for Dm, set starting values for P and W
  if (mE > 0) {
    eigenD <- eigen(Matrix::solve(Dm), symmetric = TRUE)
    Ue <- as(eigenD$vectors[, 1:mE], "dgeMatrix")
    psiE <- mean(eigenD$values[-(1:mE)])
    if (mE > 1) {
      rootLambdaE <- matrixRoot(
        Matrix::Diagonal(x = eigenD$values[1:mE] - psiE))
    } else {
      rootLambdaE <- as(sqrt(eigenD$values[1:mE] - psiE), "dgeMatrix")
    }
    We <- Ue %*% rootLambdaE
    Pe <- Matrix::Diagonal(n = p, x = 1 / psiE)
  } else {
    We <- NULL
    Pe <- NULL
  }
  ## Set starting values.
  continue <- TRUE
  decreased <- FALSE
  iter <- 1
  ELogLikCm <- ELogLikDm <- -Inf
  mu <- Matrix::Matrix(data = 0, nrow = n, ncol = p)
  ## EM following the notation of Dahl et al.
  while (continue && iter < maxIter) {
    DmSqrtInv <- matrixRoot(Matrix::solve(Dm))
    w1 <- eigen(DmSqrtInv %*% Cm %*% DmSqrtInv, symmetric = TRUE)
    Q1 <- as(w1$vectors, "dgeMatrix")
    lambda1 <- w1$values
    CmSqrtInv <- matrixRoot(Matrix::solve(Cm))
    w2 <- eigen(CmSqrtInv %*% Dm %*% CmSqrtInv, symmetric = TRUE)
    Q2 <- as(w2$vectors, "dgeMatrix")
    lambda2 <- w2$values
    if (nc > 0) {
      tUYminXb <- Matrix::crossprod(Uk, Y - X %*% B)
      S1 <- vecInvDiag(x = lambda1, y = w$values) *
        (tUYminXb %*% matrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) *
        (tUYminXb %*% matrixRoot(Cm) %*% Q2)
    } else {
      S1 <- vecInvDiag(x = lambda1, y = w$values) *
        Matrix::crossprod(Uk, Y %*% matrixRoot(Dm) %*% Q1)
      S2 <- vecInvDiag(x = lambda2, y = 1 / w$values) *
        Matrix::crossprod(Uk, Y %*% matrixRoot(Cm) %*% Q2)
    }
    trP1 <- tracePInvDiag(x = lambda1, y = w$values)
    trP2 <- tracePInvDiag(x = lambda2, y = 1 / w$values)
    if (p > 1) {
      part1 <- DmSqrtInv %*% Q1 %*% Matrix::Diagonal(x = trP1) %*%
        Matrix::crossprod(Q1, DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% Matrix::Diagonal(x = trP2) %*%
        Matrix::crossprod(Q2, CmSqrtInv)
    } else {
      part1 <- DmSqrtInv %*% Q1 %*% Matrix::Matrix(trP1) %*%
        Matrix::crossprod(Q1, DmSqrtInv)
      part2 <- CmSqrtInv %*% Q2 %*% Matrix::Matrix(trP2) %*%
        Matrix::crossprod(Q2, CmSqrtInv)
    }
    part3 <- Matrix::tcrossprod(CmSqrtInv %*% Matrix::tcrossprod(Q2, S2))
    part4 <- Matrix::tcrossprod(DmSqrtInv %*% Q1, S1) %*% lambdaR %*%
      Matrix::tcrossprod(S1, DmSqrtInv %*% Q1)
    if (nc > 0) {
      mu <- Matrix::tcrossprod(Uk %*% S1, DmSqrtInv %*% Q1)
      B <- XtXinvXt %*% (Y - mu)
    }
    Omega1 <- Matrix::forceSymmetric((part1 + part3) / n)
    Omega2 <- Matrix::forceSymmetric((part2 + part4) / n)
    ## Update Cm
    if (mG == 0) {
      ## Recall that the model is Cm^{-1} = P^{-1} + W W^t.
      ## when mG == 0, W = 0 and Cm = P
      if (p > 1) {
        if (CmHet) {
          PgNew <- Matrix::Diagonal(x = pmin(maxDiag, 1 / Matrix::diag(Omega2)))
        } else {
          tau <- min(maxDiag, p / sum(Matrix::diag(Omega2)))
          PgNew <- Matrix::Diagonal(n = p, x = tau)
        }
      } else {
        PgNew <- Matrix::Matrix(1 / as.numeric(Omega2))
      }
      WgNew <- NULL
      CmNew  <- PgNew
    } else {
      ## When rank(Omega) = Q, A should be the Q x p matrix such that
      ## Omega = A^t A / Q
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
      CmNew <- Matrix::solve(Matrix::solve(PgNew) + Matrix::tcrossprod(WgNew))
    }
    ## Update Dm
    if (mE == 0) {
      ## Recall that the model is Dm^{-1} = P^{-1} + W W^t.
      ## when mE == 0, W = 0 and Cm = P
      if (p > 1) {
        if (DmHet) {
          PeNew <- Matrix::Diagonal(x = pmin(maxDiag, 1 / diag(Omega1)))
        } else {
          tau <- min(maxDiag, p / sum(Matrix::diag(Omega1)))
          PeNew <- Matrix::Diagonal(n = p, x = tau)
        }
      } else {
        PeNew <- Matrix::Matrix(1 / as.numeric(Omega1))
      }
      WeNew <- NULL
      DmNew  <- PeNew
    } else {
      ## When rank(Omega) = Q, A should be the Q x p matrix such that
      ## Omega = A^t A / Q
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
      DmNew <- Matrix::solve(Matrix::solve(PeNew) + Matrix::tcrossprod(WeNew))
    }
    ## Compute log-likelihood and check stopping criteria
    ELogLikOld <- ELogLikCm + ELogLikDm
    ELogLikCm <- n * (Matrix::determinant(Cm)[[1]][1] -
                        sum(Matrix::diag(Cm %*% Omega2)))
    ELogLikDm <- n * (Matrix::determinant(Dm)[[1]][1] -
                        sum(Matrix::diag(Dm %*% Omega1)))
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
    VInvLst <- makeVInvLst(Vg <- Matrix::solve(Cm), Ve <- Matrix::solve(Dm),
                             Dk = Dk)
    VLst <- makeVLst(Vg = Vg, Ve = Ve, Dk = Dk)
    if (nc > 0) {
      XTransformed <- Matrix::crossprod(X, Uk)
    } else {
      XTransformed <- data.frame()
    }
    logLik <- LLDiag(Matrix::crossprod(Y, Uk), X = XTransformed,
                     vLst = VLst, vInvLst = VInvLst)
  } else {
    logLik <- NA
  }
  ## Add default names if needed.
  if (is.null(rownames(Y))) {
    rownames(Y) <- paste0("genotype", 1:n)
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("trait", 1:p)
  }
  if (prediction & nc == 0) {
    mu <- Matrix::tcrossprod(Uk %*% S1, DmSqrtInv %*% Q1)
  }
  predFrame <- data.frame(trait = rep(colnames(Y), each = n),
                          genotype = rep(rownames(Y), p),
                          predicted = as.numeric(mu))
  Vg = Matrix::solve(Cm)
  Ve = Matrix::solve(Dm)
  colnames(Vg) <- rownames(Vg) <- colnames(Ve) <- rownames(Ve) <- colnames(Y)
  return(list(Vg = Vg, Ve = Ve))
}

#' Update W and P in EMFA algorithm
#'
#' Update W and P used in the iteration process in the EMFA algorithm.
#'
#' @param Y An n x p matrix or data.frame.
#' @param WStart A p x p matrix or data.frame containing starting values for W.
#' @param m An integer. The order of the model.
#' @param PStart A p x p matrix or data.frame containing starting values for P.
#' @param hetVar Should an extra diagonal part be added in the model for the
#' precision matrix?
#' @param maxDiag A numerical value for the maximum value of the diagonal of P.
#' @param tolerance A numerical value. The iterating process stops if the sum
#' of the difference for P and W between two steps gets lower than this value.
#' @param maxIter A numerical value for the maximum number of iterations.
#' @param printProgress Should progress be printed during iterations?
#'
#' @return A list containing the new matrices W and P after the iteration
#' process and the number of iterations.
#'
#' @keywords internal
updateFA <- function(Y,
                     WStart = NULL,
                     m = ifelse(is.null(WStart), 2, ncol(WStart)),
                     PStart = NULL,
                     hetVar = FALSE,
                     maxDiag = 1e4,
                     tolerance = 1e-4,
                     maxIter = 100L,
                     printProgress = FALSE) {
  ## Check input
  if (anyNA(Y)) {
    stop("Y cannot contain missing values.\n")
  }
  p <- ncol(Y)
  n <- nrow(Y)
  if (!is.null(PStart)) {
    stopifnot(ncol(PStart) == p & nrow(PStart) == p)
  }
  if (!is.null(WStart)) {
    stopifnot(nrow(WStart) == p)
    if (ncol(WStart) != m) {
      stop("m needs to be equal to the number of columns of WStart.")
    }
    if (is.null(PStart)) {
      stop(paste("WStart and PStart should be either both NULL (default),",
                 "or both have a sensible value."))
    }
  } else {
    if (!is.null(PStart)) {
      stop(paste("WStart and PStart should be either both NULL (default),",
                 "or both have a sensible value."))
    }
  }
  if (m != round(m) || m < 1) {
    stop("m needs to be a positive integer.")
  }
  if (m >= p) {
    stop("m needs to be smaller than the number of variables.")
  }
  ## Set start values for P and W.
  if (is.null(WStart)) {
    a <- eigen(Matrix::crossprod(Y) / n, symmetric = TRUE)
    sigma2 <- mean(a$values[-(1:m)])
    PStart <- Matrix::Diagonal(n = p, x = 1 / sigma2)
    WStart <- a$vectors[, 1:m] %*%
      Matrix::Diagonal(x = sqrt(a$values[1:m] - sigma2))
  }
  W <- WStart
  P <- PStart
  ## Set start values for iterations and difference.
  iter <- 1
  totalDiff <- Inf
  ## EM
  while (totalDiff > tolerance & iter < maxIter) {
    if (m == 1) {
      if (hetVar) {
        B <- as.numeric(Matrix::crossprod(W, P %*% W)) # m x m
        Sigma <- 1 / (1 + B) # m x m
        M1 <- Sigma *
          as.numeric(Matrix::crossprod(W, Matrix::tcrossprod(P, Y))) # m x n
      } else {
        B <- P[1, 1] * as.numeric(Matrix::crossprod(W)) # m x m
        Sigma <- 1 / (1 + B) # m x m
        M1 <- P[1, 1] * Sigma * as.numeric(Matrix::t(Y %*% W)) # m x n
      }
      A <- 1 / (n * Sigma + sum(M1 ^ 2))  # m x m
      WNew <- A * Matrix::crossprod(Y, M1) # p x m
      if (hetVar) {
        D1 <- Matrix::colSums(Y ^ 2)
        D2 <- (n * Sigma + sum(M1 ^ 2)) * as.numeric(WNew) ^ 2
        D3 <- Matrix::diag(WNew %*% M1 %*% Y)
        DTot <-  D1 + D2 - 2 * D3
        PNew <- P
        Matrix::diag(PNew) <- n / DTot
      } else {
        dataNew <- Matrix::t(Y) - WNew %*% M1 # p x n
        SNew <- Matrix::tcrossprod(dataNew) / n
        PNew <- P
        Matrix::diag(PNew) <- 1 / mean(Matrix::diag(SNew))
      }
    } else {
      if (hetVar) {
        B <- Matrix::crossprod(W, P %*% W) # m x m
        Sigma <- Matrix::solve(Matrix::Diagonal(n = m) + B) # m x m
        M1 <- Sigma %*% Matrix::crossprod(W, Matrix::tcrossprod(P, Y)) # m x n
      } else {
        B <- P[1, 1] * Matrix::crossprod(W) # m x m
        Sigma <- Matrix::solve(Matrix::Diagonal(n = m) + B) # m x m
        M1 <- P[1, 1] *
          Matrix::tcrossprod(Matrix::tcrossprod(Sigma, W), Y) # m x n
      }
      A <- Matrix::solve(n * Sigma + Matrix::tcrossprod(M1))  # m x m
      WNew <- Matrix::crossprod(Y, Matrix::crossprod(M1, A)) # p x m
      if (hetVar) {
        D1 <- Matrix::colSums(Y ^ 2)
        D2 <- Matrix::diag(
          Matrix::tcrossprod(WNew %*% (n * Sigma + Matrix::tcrossprod(M1)),
                             WNew))
        D3 <- Matrix::diag(WNew %*% M1 %*% Y)
        DTot <-  D1 + D2 - 2 * D3
        PNew <- P
        Matrix::diag(PNew) <- n / DTot
      } else {
        dataNew <- Matrix::t(Y) - WNew %*% M1 # p x n
        SNew <- Matrix::tcrossprod(dataNew) / n
        PNew <- P
        Matrix::diag(PNew) <- 1 / mean(Matrix::diag(SNew))
      }
    }
    Matrix::diag(PNew)[Matrix::diag(PNew) > maxDiag] <- maxDiag
    PDiff <- sum(abs(as.numeric(PNew) - as.numeric(P)))
    WDiff <- sum(abs(as.numeric(WNew) - as.numeric(W)))
    totalDiff <- PDiff + WDiff
    if (printProgress) {
      cat("Iteration ", iter, " : ", PDiff, "  ", WDiff, "\n")
    }
    ## Set values for next iteration
    P <- PNew
    W <- WNew
    iter <- iter + 1
  }
  return(list(W = W, P = P, n.iter = iter))
}

#' Update W and P in EMFA algorithm for homogeneous variance.
#'
#' Updata W and P used in the iteration process in the EMFA algorithm in case
#' the variance is homogeneous.
#'
#' @inheritParams EMFA
#'
#' @param S A p x p sample covariance matrix.
#' @param m An integer. The order of the model.
#' @param maxDiag A numerical value for the maximum value of sigma2.
#'
#' @return A list containing the new value for the matrices W and P.
#'
#' @keywords internal
updateFAHomVar <- function(Y = NULL,
                           S = NULL,
                           m,
                           maxDiag = 1e4) {
  if ((is.null(Y) && is.null(S)) || (!is.null(Y) && !is.null(S))) {
    stop(paste("Either the data (Y) or the sample covariance matrix (S)",
               "must be provided.\n"))
  }
  if (m != round(m) || m < 1) {
    stop("m needs to be a positive integer")
  }
  ## If S is not in input, compute S from Y.
  if (!is.null(Y)) {
    if (anyNA(Y)) {
      stop("Y cannot contain missing values.\n")
    }
    n <- nrow(Y)
    Y <- as(scale(Y, scale = FALSE), "dgeMatrix")
    S <- Matrix::crossprod(Y) / n
  }
  p <- ncol(S)
  a <- eigen(S, symmetric = TRUE)
  if (m >= p) {
    stop("m needs to be smaller than the number of variables.\n")
  }
  sigma2 <- max(mean(a$values[-(1:m)]), 1 / maxDiag)
  ## Split cases for extra robustness.
  if (m == 1) {
    W <- as(a$vectors[, 1] * sqrt(a$values[1] - sigma2), "dgeMatrix")
  } else {
    W <- a$vectors[, 1:m] %*% Matrix::Diagonal(x = sqrt(a$values[1:m] - sigma2))
  }
  return(list(W = W, P = Matrix::Diagonal(n = p, x = 1 / sigma2)))
}

#' Compute log-likelihood
#'
#' Compute \eqn{t(y) * P * y}, part of the log-likelihood functions from
#' equation 26 and 27 in Zhou and Stephens (2014) using equation 50.
#' Equation 56, 57 and 58 are used to do the actual computations.
#'
#' It is assumed that X and Y have already been rotated by Uk, where Uk is such
#' that the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
#' The original X and Y are right multiplied by Uk, e.g. \code{Y <- Y * Uk}.
#' See Zhou and Stephens (2014), supplement.\cr
#' It is these rotated versions that are the input of this function.
#'
#' @inheritParams estEffs
#'
#' @param X An optional c x n covariate matrix, c being the number of
#' covariates and n being the number of genotypes. c has to be at least one
#' (typically an intercept). No missing values are allowed.
#' @param vLst A list of n pxp matrices obtained as an output from the
#' function \code{\link{makeVLst}}. It contains for each genotype l the
#' p x p matrix \eqn{v_l} (in the notation of Zhou and Stephens)
#' @param vInvLst A list of n pxp matrices obtained as an output from the
#' function \code{\link{makeVInvLst}}. It contains for each genotype l the
#' p x p matrix \eqn{v_l ^ {-1}} (in the notation of Zhou and Stephens).
#'
#' @return A numerical value for the \eqn{t(y) * P * y} part of the
#' log-likelihood function.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal
LLDiag <- function(Y,
                   X = data.frame(),
                   vLst,
                   vInvLst) {
  nc <- nrow(X)
  n <- ncol(Y)
  p <- nrow(Y)
  ## Define function for faster computation of scalar part.
  scalFunc <- function(i) {
    as.numeric(Matrix::crossprod(Y[, i, drop = FALSE],
                                 vInvLst[[i]] %*% Y[, i, drop = FALSE]))
  }
  ## Compute scalair part.
  qScal <- sum(sapply(X = 1:n, FUN = scalFunc))
  quadFormPart <- -0.5 * qScal
  if (nc > 0) {
    ## Define functions for faster computation of q and Q.
    qVecFunc <- function(i) {
      as.numeric(Matrix::kronecker(X[, i, drop = FALSE],
                                   vInvLst[[i]] %*% Y[, i, drop = FALSE]))
    }
    qMatFunc <- function(i) {
      as.numeric(Matrix::kronecker(Matrix::tcrossprod(X[, i, drop = FALSE]),
                                   vInvLst[[i]]))
    }
    ## Compute q, Q and quadratic part.
    qVec <- rowSums(sapply(X = 1:n, FUN = qVecFunc))
    QMatrix <- matrix(rowSums(sapply(X = 1:n, FUN = qMatFunc)), ncol = p * nc)
    quadFormPart <- quadFormPart + 0.5 *
      as.numeric(crossprod(qVec, solve(QMatrix, qVec)))
  }
  ## Compute determinant part.
  detPart <- -0.5 * sum(sapply(X = 1:n, FUN = function(i) {
    Matrix::determinant(vLst[[i]])[[1]][1]
  }))
  return(quadFormPart + detPart)
}

#' Compute tYPY as in Zhou and Stephens eqn. 50.
#'
#' Compute \eqn{t(y) * P * y}, part of the log-likelihood functions from
#' equation 26 and 27 in Zhou and Stephens using equation 50. Equation 56, 57
#' and 58 are used to do the actual computations.
#'
#' It is assumed that X and Y have already been rotated by Uk, where Uk is such
#' that the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
#' The original X and Y are right multiplied by Uk, e.g. \code{Y <- Y * Uk}.
#' See Zhou and Stephens (2014), supplement.\cr
#' It is these rotated versions that are the input of this function.
#'
#' @inheritParams LLDiag
#'
#' @return A numerical value for the \eqn{t(y) * P * y} part of the
#' log-likelihood function.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal
LLQuadFormDiag <- function(Y,
                           X = data.frame(),
                           vInvLst) {
  nc <- nrow(X)
  n <- ncol(Y)
  p <- nrow(Y)
  ## Define function for faster computation of scalar part.
  scalFunc <- function(i) {
    as.numeric(Matrix::crossprod(Y[, i, drop = FALSE],
                                 vInvLst[[i]] %*% Y[, i, drop = FALSE]))
  }
  ## Compute scalair part.
  qScal <- sum(sapply(X = 1:n, FUN = scalFunc))
  if (nc > 0) {
    ## Define functions for faster computation of q and Q.
    qVecFunc <- function(i) {
      as.numeric(Matrix::kronecker(X[, i, drop = FALSE],
                                   vInvLst[[i]] %*% Y[, i, drop = FALSE]))
    }
    qMatFunc <- function(i) {
      as.numeric(Matrix::kronecker(Matrix::tcrossprod(X[, i, drop = FALSE]),
                                   vInvLst[[i]]))
    }
    ## Compute q and Q.
    if (p == 1 && nc == 1) {
      qVec <- sum(sapply(X = 1:n, FUN = qVecFunc))
      QMatrix <- sum(sapply(X = 1:n, FUN = qMatFunc))
    } else {
      qVec  <- rowSums(sapply(X = 1:n, FUN = qVecFunc))
      QMatrix <- matrix(rowSums(sapply(X = 1:n, FUN = qMatFunc)), ncol = p * nc)
    }
    ## Compute quadratic part.
    quadFormPart <- qScal - as.numeric(crossprod(qVec %*% solve(QMatrix, qVec)))
  } else {
    quadFormPart <- qScal
  }
  return(quadFormPart)
}

#' Create a list of (inverted) variance matrices
#'
#' As in Zhou and Stephens, supplement page 13, create the list of (inverted)
#' variance matrices for the transformed genotypes per individual as given by
#' equation 4.
#'
#' @param Vg A p x p symmetric matrix of genetic variance components
#' @param Ve A p x p symmetric matrix of environmental variance components
#' @param Dk A vector of length n containing the eigenvalues obtained by the
#' eigen-decomposition of the kinship matrix K: \eqn{K = Uk * Dk * t(Uk)}
#'
#' @return A list of n p x p matrices \eqn{v_l} where \eqn{v_l =
#' Dk_l * Vg + Ve \forall l = 1, ..., n}.\cr
#' When using \code{makeVInvLst} the output matrices are inverted.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409.
#'
#' @keywords internal
makeVLst <- function(Vg,
                     Ve,
                     Dk) {
  VLst <- sapply(X = Dk, FUN = function(Dk_i) {
    Dk_i * Vg + Ve
  }, simplify = FALSE)
  return(VLst)
}

#' @rdname makeVLst
makeVInvLst <- function(Vg,
                        Ve,
                        Dk) {
  VInvLst <- sapply(X = Dk, FUN = function(Dk_i) {
    Matrix::solve(Dk_i * Vg + Ve)
  }, simplify = FALSE)
  return(VInvLst)
}

#' Helper functions for the penalized EM algorithm
#'
#' \code{vecInvDiag} is a helper function for quickly computing
#' \eqn{(I + x \otimes y)^{-1}},
#' \code{tracePInvDiag} for quickly computing column sums of
#' \eqn{(I + x \otimes y)^{-1}}. Both are used in the penalized EM algorithm.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#'
#' @return for \code{vecInvDiag} a matrix defined by
#' \eqn{(I + x \otimes y)^{-1}}, for \code{tracePInvDiag} a vector containing
#' the column sums of \eqn{(I + x \otimes y)^{-1}}.
#'
#' @keywords internal
vecInvDiag <- function(x, y) {
  z <- sapply(X = x, FUN = function(x_i) {
    1 / (1 + x_i * y)
  })
  return(z)
}

#' @rdname vecInvDiag
tracePInvDiag <- function(x, y) {
  z <- sapply(X = x, FUN = function(x_i) {
    sum(1 / (1 + x_i * y))
  })
  return(z)
}
