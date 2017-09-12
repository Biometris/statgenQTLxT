#' Update W and P in EMFA algorithm
#'
#' Updata W and P used in the iteration process in the EMFA algorithm.
#'
#' @param Y a n x p matrix or dataframe.
#' @param WStart a p x p matrix or dataframe containing starting values for W.
#' @param m an integer. The order of the model.
#' @param PStart a p x p matrix or dataframe containing starting values for P.
#' @param hetVar should an extra diagonal part be added in the model for the precision matrix?
#' @param maxDiag a numerical value for the maximum value of the diagonal of P.
#' @param tolerance a numerical value. The iterating process stops if the sum of the difference for P
#' and W between two steps gets lower than this value.
#' @param maxIter a numerical value for the maximum number of iterations.
#' @param printProgress should progress be printed during iterations?
#'
#' @return a list containing the new matrices W and P after the iteration process and the number
#' of iterations.
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
  Y <- as.matrix(Y)
  if (anyNA(Y)) {stop('Y cannot contain missing values')}
  p <- ncol(Y)
  n <- nrow(Y)
  if (!is.null(PStart)) {
    stopifnot(class(PStart) %in% c('matrix','data.frame'))
    stopifnot(ncol(PStart) == p & nrow(PStart) == p)
  }
  if (!is.null(WStart)) {
    stopifnot(class(WStart) %in% c('matrix','data.frame'))
    stopifnot(nrow(WStart) == p)
    if (ncol(WStart) != m) {stop('m needs to be equal to the number of columns of WStart')}
    if (is.null(PStart)) {stop('WStart and PStart should be either both NULL (default),
      or both have a sensible value')}
  } else {
    if (!is.null(PStart)) {stop('WStart and PStart should be either both NULL (default),
      or both have a sensible value')}
  }
  if (m != round(m) || m < 1) {stop("m needs to be a positive integer")}
  if (m >= p) {stop("m needs to be smaller than the number of variables")}
  ## Set start values for P and W.
  if (is.null(WStart)) {
    a <- eigen(crossprod(Y) / n, symmetric = TRUE)
    sigma2 <- mean(a$values[-(1:m)])
    PStart <- diag(x = 1 / sigma2, nrow = p)
    WStart <- a$vectors[, 1:m] %*% diag(sqrt(a$values[1:m] - sigma2))
  }
  W <- WStart
  P <- PStart
  ## Set start values for iterations and difference.
  iter <- 1
  totalDiff <- Inf
  ## EM
  while (totalDiff > tolerance & iter < maxIter) {
    ## prevent that P become asymmetric because of numerical inaccuracies
    P <- (P + t(P)) / 2 # p x p
    if (hetVar) {
      B <- crossprod(W, P %*% W) # m x m
      sigma <- MASS::ginv(diag(m) + B) # m x m
      M1 <- sigma %*% crossprod(W, tcrossprod(P, Y)) # m x n
    } else {
      B <- P[1, 1] * crossprod(W) # m x m
      sigma <- MASS::ginv(diag(m) + B) # m x m
      M1 <- P[1, 1] * tcrossprod(tcrossprod(sigma, W), Y) # m x n
    }
    A <- MASS::ginv(n * sigma + tcrossprod(M1))  # m x m
    WNew <- crossprod(Y, crossprod(M1, A)) # p x m
    if (hetVar) {
      D1 <- colSums(Y ^ 2)
      D2 <- diag(WNew %*% (n * sigma + tcrossprod(M1)) %*% t(WNew))
      D3 <- diag(WNew %*% M1 %*% Y)
      DTot <-  D1 + D2 - 2 * D3
      PNew <- diag(n / DTot)
    } else {
      dataNew <- t(Y) - WNew %*% M1 # p x n
      SNew <- tcrossprod(dataNew) / n
      PNew <- diag(1 / mean(diag(SNew)), nrow = ncol(SNew))
    }
    diag(PNew)[diag(PNew) > maxDiag] <- maxDiag
    PDiff <- sum(abs(PNew - P))
    WDiff <- sum(abs(WNew - W))
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

