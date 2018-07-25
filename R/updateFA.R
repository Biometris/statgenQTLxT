#' Update W and P in EMFA algorithm
#'
#' Updata W and P used in the iteration process in the EMFA algorithm.
#'
#' @param Y a n x p matrix or dataframe.
#' @param WStart a p x p matrix or dataframe containing starting values for W.
#' @param m an integer. The order of the model.
#' @param PStart a p x p matrix or dataframe containing starting values for P.
#' @param hetVar should an extra diagonal part be added in the model for the
#' precision matrix?
#' @param maxDiag a numerical value for the maximum value of the diagonal of P.
#' @param tolerance a numerical value. The iterating process stops if the sum
#' of the difference for P and W between two steps gets lower than this value.
#' @param maxIter a numerical value for the maximum number of iterations.
#' @param printProgress should progress be printed during iterations?
#'
#' @return a list containing the new matrices W and P after the iteration
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

