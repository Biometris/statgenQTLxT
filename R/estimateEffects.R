#' Estimates for covariates.
#'
#' Compute the estimates and standard errors for the covariates in the input matrix X. Optionally an
#' additional covariate x can be estimated.
#'
#' No missing values are allowed in X, Y and x.\cr
#' It is assumed that X, Y and x have already been rotated by Uk, where Uk is such that
#' the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
#' The original X, Y and x are right multiplied by Uk, e.g. \code{Y <- Y * Uk}. See Zhou
#' and Stephens 2014, supplement.\cr
#' It is these rotated versions that are the input of this function.
#'
#' @export
#'
#' @param X a c x n covariate matrix, c being the number of covariates and n being the number
#' of genotypes. c has to be at least one (typically an intercept). See details.
#' @param Y a p x n matrix of observed phenotypes, on p traits or environments for n genotypes.
#' See details.
#' @param x an optional additional covariate, a numeric vector of length n. See details.
#' @param VInvArray an n x p x p dimensional array obtained as an output from the function
#' \code{\link{makeVInvArray}}. It contains for each genotype l the p x p matrix \eqn{v_l ^ {-1}} (in
#' the notation of Zhou and Stephens)
#' @param returnAllEffects If \code{FALSE} only the p effect estimates for the additional covariate
#' are returned, (the first p * c are actually computed, but not returned). If \code{TRUE},
#' also the effect estimates and standard errors for the other covariates (i.e. those in X) are returned.
#' If no x is provided, \code{returnAllEffects} is automatically set to \code{TRUE}.
#'
#' @return A list containing the estimates, the standard errors of the estimates and if an value for
#' x is input
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear mixed model algorithms for
#' genome-wide association studies.

## to do: example
## comment on missing values
## partial update...
## wald?

estimateEffects <- function(X,
  x = NULL,
  Y,
  VInvArray,
  returnAllEffects = TRUE) {

  stopifnot(ncol(X) == ncol(Y))
  stopifnot(ncol(X) == dim(VInvArray)[1])
  stopifnot(nrow(Y) == dim(VInvArray)[2])
  stopifnot(nrow(Y) == dim(VInvArray)[3])
  stopifnot(is.null(x) || length(x) == ncol(X))
  if (anyNA(X)) stop("No missing values allowed in X")
  if (anyNA(x)) stop("No missing values allowed in x")
  if (anyNA(Y)) stop("No missing values allowed in Y")

  nc <- nrow(X)
  n <- ncol(X)
  p <- nrow(Y)

  if (is.null(x)) {
    returnAllEffects <- TRUE
    ncTot <- nc
  } else {
    ncTot <- nc + 1
  }

  ## the last p coefficients should be the marker-effects
  ## the first p * nc should correspond to the other coefficients
  if (!is.null(x)) {X <- rbind(X , t(matrix(as.numeric(x))))}

  if (p == 1 && ncTot == 1) {
    Vbeta <- matrix(sum(sapply(1:n, function(i) {
      kronecker(matrix(X[, i]) %*% t(matrix(X[, i])), VInvArray[i, , ])})))
  } else {
    Vbeta <- matrix(rowSums(sapply(1:n, function(i) {
      kronecker(matrix(X[, i]) %*% t(matrix(X[, i])), VInvArray[i, , ])})),
      ncol = p * ncTot)
  }

  M <- solve(Vbeta)

  if (length(x) > 0) {
    MSub <- solve(as.matrix(Vbeta[-(1:(p * nc)), -(1:(p * nc))]))
  }

  if (p == 1 && ncTot == 1) {
    v <- sum(sapply(1:n, function(i) {kronecker(matrix(X[, i]), VInvArray[i, , ] %*% matrix(Y[, i]))}))
  } else {
    v <- rowSums(sapply(1:n, function(i) {kronecker(matrix(X[, i]), VInvArray[i, , ] %*% matrix(Y[, i]))}))
  }

  wald  <- NA

  if (returnAllEffects) {
    effectsEstimates <- as.numeric(M %*% v)
    effectsSd <- sqrt(diag(M))
    if (length(x) > 0) {
      wald <- as.numeric(t(matrix(effectsEstimates[-(1:(p * nc))])) %*% MSub %*%
          matrix(effectsEstimates[-(1:(p * nc))]))
    }
  } else {
    effectsEstimates <- as.numeric(M %*% v)[-(1:(p * nc))]
    effectsSd <- sqrt(diag(M))[-(1:(p * nc))]
    if (length(x) > 0) {
      wald <- as.numeric(t(matrix(effectsEstimates)) %*% MSub %*% matrix(effectsEstimates))
    }
  }
  return(list(effects.estimates = effectsEstimates, effects.sd = effectsSd, wald = wald))
}
