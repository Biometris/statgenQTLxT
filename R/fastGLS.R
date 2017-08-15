#' fastGLS
#'
#' compute p-values for the GLS F-test as in emma-x. \code{fastGLS} can be used when
#' there are no further covariates, otherwise use \code{fastGLSCov}.
#'
#' @param y a numeric vector of length n of phenotypic scores. No missing values allowed.
#' @param X an n x m matrix of marker-scores, n being the number of individuals, m the number of markers.
#' no missing values allowed.
#' @param Sigma an n x n covariance matrix. No missing values allowed.
#' @param covs an n x c matrix of covariates (NOT including an intercept). No missing values allowed.
#' @param nChunks an integer, the number of parts in which the calculations are split.
#'
#' @return a data.frame with the following columns:
#' \itemize{
#' \item{\code{pValue} a vector of p-values for the GLS F-test}
#' \item{\code{beta} a vector of effect sizes}
#' \item{\code{betaSe} a vector of standard errors of the effect sizes}
#' \item{\code{RLR2} a vector of R_LR^2 statistics as defined in Sun et al.}
#' }
#'
#' @references Segura et al. (2012) An efficient multi-locus mixed-model approach for genome-wide
#' association studies in structured populations. Nature Genetics, June 2012, Vol. 44, p. 825–830.
#' @references Sun et al. (2010) Variation explained in mixed-model association mapping.
#' Heredity, February 2010, Vol. 105, p. 333–340.
#'
#' @import stats

fastGLS <- function(y,
  X,
  Sigma) {
  ## Check class and missing values.
  if (missing(y) || !is.numeric(y) || anyNA(y))
    stop("y should be a numeric vector without missing values.")
  if (missing(X) || !is.matrix(X) || !is.numeric(X) || anyNA(X))
    stop("y should be a matrix without missing values.")
  if (missing(Sigma) || !is.matrix(Sigma) || !is.numeric(Sigma) || anyNA(Sigma))
    stop("Sigma should be a matrix without missing values.")
  ## Check dimensions.
  n <- length(y)
  if (nrow(X) != n)
    stop("The number of elements in y should be identical to the number of rows in X")
  if (nrow(Sigma) != n || ncol(Sigma) != n)
    stop("The number of elements in y should be identical to the number of rows and columns in Sigma")
  # check for missing values, and class of y
  M <- solve(chol(Sigma))
  ## pre-multiply the phenotype (y) with t(M)
  tMy <- crossprod(M, y)
  ## pre-multiply the intercept with t(M)
  tMInt <- crossprod(M, rep(1, n))
  ## pre-multiply the snp-matrix with t(M)
  ## for extra robustness, distinguish
  if (ncol(X) == 1) {
    tMX <- crossprod(M, matrix(as.numeric(X)))
  } else {
    tMX <- crossprod(M, X)
  }
  ## Matrix cookbook, 3.2.6 Rank-1 update of inverse of inner product
  a  <- 1 / as.numeric(crossprod(tMInt))
  vv <- colSums(tMX ^ 2)
  vX <- as.numeric(crossprod(tMInt,  tMX))
  nn <- 1 / (vv - a * vX ^ 2)
  XtXinv2ndRows <- cbind(-1 * a * vX * nn, nn)
  Xty <- cbind(rep(as.numeric(crossprod(tMInt, tMy), length(nn))), as.numeric(crossprod(tMy, tMX)))
  betaVec <- XtXinv2ndRows[, 1] * Xty[, 1] + XtXinv2ndRows[, 2] * Xty[, 2]
  ## Compute residuals and RSS over all markers.
  RSSEnv <- sum(lsfit(x = tMInt, y = tMy, intercept = FALSE)$residuals ^ 2)
  ## Compute RSS per marker.
  RSSFull <- apply(tMX, 2, function(x) {
    sum(lsfit(x = cbind(tMInt, x), y = tMy, intercept = FALSE)$residuals ^ 2)})
  ## Compute F and p values.
  df2 <- n - 2
  FVal <- (RSSEnv - RSSFull) / RSSFull * df2
  pVal <- pf(q = FVal, df1 = 1, df2 = df2, lower.tail = FALSE)
  # the R_LR^2 statistic from G. Sun et al 2010, heredity
  RLR2  <- 1 - exp((RSSFull - RSSEnv) / n)
  ## Construct output data.frame.
  GLS <- data.frame(pValue = pVal,
    beta = betaVec,
    betaSe = sqrt(XtXinv2ndRows[, 2]),
    RLR2 = RLR2)
  return(GLS)
}

#' @rdname fastGLS
fastGLSCov <-function(y,
  X,
  Sigma,
  covs,
  nChunks = 10) {
  ## Check class and missing values.
  if (missing(y) || !is.numeric(y) || anyNA(y))
    stop("y should be a numeric vector without missing values.")
  if (missing(X) || !is.matrix(X) || !is.numeric(X) || anyNA(X))
    stop("y should be a matrix without missing values.")
  if (missing(Sigma) || !is.matrix(Sigma) || !is.numeric(Sigma) || anyNA(Sigma))
    stop("Sigma should be a matrix without missing values.")
  if (missing(covs) || !is.numeric(covs) || anyNA(covs))
    stop("covs should be a numeric vector without missing values.")
  if (!is.numeric(nChunks) || length(nChunks) > 1 || nChunks != round(nChunks))
    stop("nChunks should be an integer")
  n <- length(y)
  ## Check dimensions.
  if (nrow(X) != n)
    stop("The number of elements in y should be identical to the number of rows in X")
  if (nrow(Sigma) != n || ncol(Sigma) != n)
    stop("The number of elements in y should be identical to the number of rows and columns in Sigma")
  if (nrow(covs) != n)
    stop("The number of elements in y should be identical to the number of rows in covs")
  m <- ncol(X)
  fixCovs <- cbind(rep(1, n), covs)
  nCov <- ncol(fixCovs)
  M <- solve(chol(Sigma))
  ## pre-multiply the phenotype (y) with t(M)
  tMy <- crossprod(M, y)
  ## pre-multiply the intercept and covariates with t(M)
  tMfixCovs <- crossprod(M, fixCovs)
  ## pre-multiply the snp-matrix with t(M)
  ## for extra robustness, distinguish
  if (m == 1) {
    tMX <- crossprod(M, matrix(as.numeric(X)))
  } else {
    tMX <- crossprod(M, X)
  }
  ## Matrix cookbook, 3.2.6 Rank-1 update of inverse of inner product.
  A <- solve(crossprod(tMfixCovs), symmetric = TRUE)
  vv <- colSums(tMX ^ 2)
  vX <- crossprod(tMfixCovs, tMX)
  nn <- 1 / (vv - colSums(vX * (A %*% vX)))
  XtXinvLastRows <- cbind(- nn * crossprod(vX, A), nn)
  Xty <- cbind(matrix(rep(as.numeric(crossprod(tMfixCovs, tMy)), length(nn)), byrow = TRUE, ncol= nCov)
    , as.numeric(crossprod(tMy, tMX)))
  betaVec <- rowSums(XtXinvLastRows[, 1:nCov] * Xty[, 1:nCov]) +
    XtXinvLastRows[, 1 + nCov] * Xty[, 1 + nCov]
  ## Compute residuals and RSS over all markers.
  ResEnv <- lsfit(x = tMfixCovs, y = tMy, intercept = FALSE)$residuals
  RSSEnv <- sum(ResEnv ^ 2)
  ## QR decomposition of covariates.
  Q <- qr.Q(qr(tMfixCovs))
  tMQtQ <- t(M %*% (diag(n) - tcrossprod(Q)))
  ## Compute RSS per marker, breaking up X for speed.
  RSSFull <- vector(mode = "list", length = nChunks)
  for (j in 1:(nChunks - 1)) {
    tX <- tMQtQ %*% X[, ((j - 1) * round(m / nChunks) + 1):(j * round(m / nChunks))]
    RSSFull[[j]] <- apply(tX, 2, function(x){
      sum(lsfit(x = x, y = ResEnv, intercept = FALSE)$residuals ^ 2)})
  }
  tX <- tMQtQ %*%  X[ , -(1:(j * round(m / nChunks)))]
  RSSFull[[nChunks]] <- apply(tX, 2, function(x){
    sum(lsfit(x = x, y = ResEnv, intercept = FALSE)$residuals ^ 2)})
  RSSFull <- unlist(RSSFull)
  ## Compute F and p values.
  df2 <- n - 1 - nCov
  FVal <- (RSSEnv - RSSFull) / RSSFull * df2
  pVal <- pf(q = FVal, df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute R_LR^2 statistic from G. Sun et al 2010, heredity.
  RLR2  <- 1 - exp((RSSFull - RSSEnv) / n)
  ## Construct output data.frame.
  GLS <- data.frame(pValue = pVal,
    beta = betaVec,
    betaSe = sqrt(XtXinvLastRows[, 1 + nCov]),
    RLR2 = RLR2)
  return(GLS)
}
