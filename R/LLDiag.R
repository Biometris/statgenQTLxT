#' LLDiag
#'
#' Compute \eqn{t(y) * P * y}, part of the log-likelihood functions from equation 26 and 27 in Zhou and
#' Stephens using equation 50. Equation 56, 57 and 58 are used to do the actual computations.
#'
#' No missing values are allowed in X and Y.\cr
#' It is assumed that X and Y have already been rotated by Uk, where Uk is such that
#' the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
#' The original X and Y are right multiplied by Uk, e.g. \code{Y <- Y * Uk}. See Zhou
#' and Stephens 2014, supplement.\cr
#' It is these rotated versions that are the input of this function.
#'
#' @inheritParams estimateEffects
#'
#' @param X an optional c x n covariate matrix, c being the number of covariates and n being the number
#' of genotypes. c has to be at least one (typically an intercept). See details.
#'
#' @return a numeric value for the \eqn{t(y) * P * y} part of the log-likelihood function.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear mixed model algorithms for
#' genome-wide association studies.

## To do: example
## To do: explanation missing X

LLDiag <- function(Y,
  X = data.frame(),
  VArray,
  VInvArray) {

  stopifnot(nrow(Y) == dim(VInvArray)[2])
  stopifnot(nrow(Y) == dim(VInvArray)[3])
  if (anyNA(Y)) stop("No missing values allowed in Y")
  if (!(is.data.frame(X) || nrow(X) == 0)) {
    stopifnot(ncol(X) == dim(VInvArray)[1])
    if (anyNA(X)) stop("No missing values allowed in X")
  }

  nc <- nrow(X)
  n <- ncol(Y)
  p <- nrow(Y)

  qScal <- sum(sapply(1:n, function(i) {t(matrix(Y[, i])) %*% VInvArray[i, , ] %*% matrix(Y[, i])}))

  quadFormPart <- -0.5 * qScal

  if (nc > 0) {
    qVec  <- matrix(rowSums(sapply(1:n, function(i)
      {kronecker(matrix(X[, i]), VInvArray[i, , ] %*% matrix(Y[, i]))})))
    QMatrix <- matrix(rowSums(sapply(1:n, function(i)
      {kronecker(matrix(X[, i]) %*% t(matrix(X[, i])), VInvArray[i, , ])})), ncol = p * nc)
    quadFormPart <- quadFormPart + 0.5 * as.numeric(t(qVec) %*% solve(QMatrix) %*% qVec)
  }

  detPart <- -0.5 * sum(sapply(1:n, function(i) {determinant(VArray[i, , ])[[1]][1]}))
  MLLogLik <- quadFormPart + detPart

  return(MLLogLik)
}
