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
#' @inheritParams estimateEffects
#'
#' @param X an optional c x n covariate matrix, c being the number of
#' covariates and n being the number of genotypes. c has to be at least one
#' (typically an intercept). No missing values are allowed.
#' @param VArray an n x p x p dimensional array obtained as an output from the
#' function \code{\link{makeVArray}}. It contains for each genotype l the
#' p x p matrix \eqn{v_l} (in the notation of Zhou and Stephens)
#'
#' @return a numeric value for the \eqn{t(y) * P * y} part of the
#' log-likelihood function.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies. Nature Methods,
#' February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal
LLDiag <- function(Y,
                   X = data.frame(),
                   VArray,
                   VInvArray) {
  nc <- nrow(X)
  n <- ncol(Y)
  p <- nrow(Y)
  ## Compute scalair part.
  qScal <- sum(sapply(X = 1:n, FUN = function(i) {
    as.numeric(Matrix::crossprod(Y[, i, drop = FALSE], VInvArray[i, , ] %*%
                                   Y[, i, drop = FALSE]))
  }))
  quadFormPart <- -0.5 * qScal
  if (nc > 0) {
    ## Compute q, Q and quadratic part.
    qVec <- rowSums(sapply(X = 1:n, FUN = function(i) {
      kronecker(X[, i], VInvArray[i, , ] %*% Y[, i])
    }))
    QMatrix <- matrix(rowSums(sapply(X = 1:n, FUN = function(i) {
      kronecker(tcrossprod(X[, i]), VInvArray[i, ,])
    })), ncol = p * nc)
    quadFormPart <- quadFormPart + 0.5 *
      as.numeric(crossprod(qVec, solve(QMatrix, qVec)))
  }
  ## Compute determinant part.
  detPart <- -0.5 * sum(sapply(X = 1:n, FUN = function(i) {
    Matrix::determinant(VArray[i, , ])[[1]][1]
  }))
  return(quadFormPart + detPart)
}
