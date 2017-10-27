#' Estimates for covariates
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
#' genome-wide association studies. Nature Methods, February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal

estimateEffects <- function(X,
                            x = NULL,
                            Y,
                            VInvArray,
                            returnAllEffects = TRUE) {
  ## Checks.
  stopifnot(ncol(X) == ncol(Y))
  #stopifnot(ncol(X) == dim(VInvArray)[1])
  #stopifnot(nrow(Y) == dim(VInvArray)[2])
  #stopifnot(nrow(Y) == dim(VInvArray)[3])
  stopifnot(is.null(x) || length(x) == ncol(X))
  if (anyNA(X)) {
    stop("No missing values allowed in X.\n")
  }
  if (anyNA(x)) {
    stop("No missing values allowed in x.\n")
  }
  if (anyNA(Y)) {
    stop("No missing values allowed in Y.\n")
  }
  nc <- nrow(X)
  n <- ncol(X)
  p <- nrow(Y)
  if (is.null(x)) {
    returnAllEffects <- TRUE
    ncTot <- nc
  } else {
    ncTot <- nc + 1
  }
  ## the last p coefficients should be the marker effects.
  ## the first p * nc should correspond to the other coefficients.
  if (!is.null(x)) {
    X <- Matrix::rbind2(X, x)
  }
  ## Define functions for faster computation of Vbeta and v.
  VbetaFunc <- function(i) {
    kronecker(tcrossprod(X[, i]), VInvArray[[i]])}
  vFunc <- function(i) {
    kronecker(X[, i], VInvArray[[i]] %*% Y[, i])}
  if (p == 1 && ncTot == 1) {
    Vbeta <- sum(sapply(X = 1:n, FUN = VbetaFunc))
  } else {
    Vbeta <- Matrix::Matrix(rowSums(sapply(X = 1:n, FUN = VbetaFunc)), ncol = p * ncTot)
  }
  M <- Matrix::solve(Vbeta)
  if (length(x) > 0) {
    MSub <- Matrix::solve(Vbeta[-(1:(p * nc)), -(1:(p * nc))])
  }
  ## Split cases for extra robustness.
  if (p == 1 && ncTot == 1) {
    v <- sum(sapply(X = 1:n, FUN = vFunc))
  } else {
    v <- rowSums(sapply(X = 1:n, FUN = vFunc))
  }
  wald <- NA
  if (returnAllEffects) {
    effectsEstimates <- M %*% v
    effectsSd <- sqrt(Matrix::diag(M))
    if (length(x) > 0) {
      wald <- as.numeric(Matrix::crossprod(effectsEstimates[-(1:(p * nc))], MSub %*%
                                             effectsEstimates[-(1:(p * nc))]))
    }
  } else {
    effectsEstimates <- M %*% v[-(1:(p * nc))]
    effectsSd <- sqrt(Matrix::diag(M)[-(1:(p * nc))])
    if (length(x) > 0) {
      wald <- as.numeric(Matrix::crossprod(effectsEstimates, MSub %*% effectsEstimates))
    }
  }
  return(list(effectsEstimates = as.numeric(effectsEstimates), effectsSd = effectsSd, wald = wald))
}
