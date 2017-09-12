#' EMMA helper functions
#'
#' Helper functions for computing REML estimates of genetic and residual variance components
#' using the EMMA algorithm.
#'
#' @inheritParams runEmma
#' @param X a q x n covariate matrix, q being the number of covariates and n being the number
#' of genotypes. q has to be at least one (typically an intercept).
#'
#' @keywords internal

emmaEigenR <- function(K,
  X) {
  n <- nrow(X)
  q <- ncol(X)
  ## Compute n-q non-zero eigenvalues of SHS as defined in eqn. 5 of Kang.
  S <- diag(n) - X %*% solve(crossprod(X), t(X))
  eig <- eigen(S %*% (K + diag(n)) %*% S, symmetric = TRUE)
  if(is.complex(eig$values))
    stop("Complex eigen values found.\n")
  return(list(values = eig$values[1:(n - q)] - 1,
    vectors = eig$vectors[, 1:(n - q)]))
}

emmaEigenRZ <- function(Z,
  K,
  X,
  complete = TRUE) {
  if (!complete) {
    vIds <- colSums(Z) > 0
    Z <- Z[, vIds]
    K <- K[vIds, vIds]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  SZ <- Z - X %*% solve(crossprod(X), crossprod(X, Z))
  eig <- eigen(K %*% tcrossprod(Z, SZ), symmetric = FALSE)
  if (is.complex(eig$values)) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)
  }
  qrX <- qr.Q(qr(X))
  return(list(values = eig$values[1:(t - q)],
    vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[, 1:(t - q)], qrX)), complete = TRUE)[, c(1:(t - q), (t + 1):n)]))
}

emmaREMLLL <- function(logDelta, lambda, etas1, n, t, etas2) {
  ## Compute the REML LL as in eqn. 7 of Kang.
  nq <- length(etas1) + n - t
  delta <- exp(logDelta)
  lDelta <- lambda + delta
  return(0.5 * (nq * (log(nq / (2 * pi)) - 1 - log(sum(etas1 ^ 2 / (lDelta)) + etas2 / delta)) -
      (sum(log(lDelta)) + (n - t) * logDelta)))
}

