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
#' @inheritParams estimateEffects
#'
#' @param X An optional c x n covariate matrix, c being the number of covariates
#' and n being the number of genotypes. c has to be at least one
#' (typically an intercept). No missing values are allowed.
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
                           vInvArr) {
  nc <- nrow(X)
  n <- ncol(Y)
  p <- nrow(Y)
  ## Define functions for faster computation of scalar part, q and Q.
  scalFunc <- function(i) {
    as.numeric(crossprod(Y[, i], vInvArr[i, ,] %*% Y[, i]))
  }
  qVecFunc <- function(i) {
    kronecker(X[, i], vInvArr[i, ,] %*% Y[, i])
  }
  qMatFunc <- function(i) {
    kronecker(tcrossprod(X[, i]), vInvArr[i, ,])
  }
  ## Compute scalair part.
  qScal <- sum(sapply(X = 1:n, FUN = scalFunc))
  if (nc > 0) {
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

