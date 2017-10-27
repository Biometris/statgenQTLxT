#' Test for marker effect.
#'
#' Test the effect of an extra covariate. This is done by testing the results of a NULL model
#' only containing the trait specific means and an alternative model containing the trait
#' specific marker effect as well.
#'
#' No missing values are allowed in X, Y and x.\cr
#' It is assumed that X, Y and x have already been rotated by Uk, where Uk is such that
#' the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
#' The original X, Y and x are right multiplied by Uk, e.g. \code{Y <- Y * Uk}. See Zhou
#' and Stephens 2014, supplement.\cr
#' It is these rotated versions that are the input of this function.
#'
#' @inheritParams estimateEffects
#'
#' @param x an additional covariate, usually a single marker, a numeric vector of length n. See details.
#' @param SS0 the log-likelihood value of a NULL model. This is the model with trait specific means
#' only and without an additional marker effect. If \code{NULL} then this value is calculated.
#'
#' @return A list containing the test statistics \code{pvalue}, \code{FStat} and \code{df}. Furthermore
#' the Log-likelihood of the null model and the alternative model, \code{SS0} and \code{SS1}. Also the
#' estimated effect, its standard error and results of the wald test for the additional covariate a returned.
#'
#' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear mixed model algorithms for
#' genome-wide association studies. Nature Methods, February 2014, Vol. 11, p. 407–409
#'
#' @keywords internal

LRTTest <- function(X,
                    x,
                    Y,
                    VInvArray,
                    SS0 = NULL) {
  nc <- nrow(X)
  n <- ncol(X)
  p <- nrow(Y)
  dfFull <- (n - (nc + 1)) * p
  ## Null model with the trait specific means only
  if (is.null(SS0)) {
    est0 <- estimateEffects(X = X, Y = Y, VInvArray = VInvArray, returnAllEffects = TRUE)
    fittedMean0 <- Matrix::Matrix(est0$effectsEstimates, nrow = p) %*% X
    SS0 <- LLQuadFormDiag(Y = Y - fittedMean0, VInvArray = VInvArray)
  }
  ## Alternative model, with additional trait specific marker effect
  est1 <- estimateEffects(X = X, x = x, Y = Y, VInvArray = VInvArray, returnAllEffects = TRUE)
  fittedMean1 <- Matrix::Matrix(est1$effectsEstimates, nrow = p) %*%
    Matrix::rbind2(X, x)
  SS1 <- LLQuadFormDiag(Y = Y - fittedMean1, VInvArray = VInvArray)
  ## Compute F-statistic and p-value using results from null and alternative model.
  FStat <- ((SS0 - SS1) / SS1) * dfFull / p
  pValue <- pf(q = FStat, df1 = p, df2 = dfFull, lower.tail = FALSE)
  return(list(pvalue = pValue, Fstat = FStat, df = p,
              SS1 = SS1, SS0 = SS0,
              effects = est1$effectsEstimates[-(1:(nc * p))],
              effectsSe = est1$effectsSd[-(1:(nc * p))],
              wald = est1$wald
  ))
}
