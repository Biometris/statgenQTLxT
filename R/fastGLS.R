#' Fast Generalized Least Squares algoritm
#'
#' Compute statistics for the Generalized Least Squares (GLS) F-test based on the algorithm proposed by
#' Segura (2012). Also the \eqn{R_LR^2} statistics as in Sun (2010) is computed.
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
#'
#' @keywords internal

fastGLS <-function(y,
  X,
  Sigma,
  covs = NULL,
  nChunks = 10) {
  ## Check class and missing values.
  if (missing(y) || !(inherits(y, "Matrix") || is.numeric(y)) || anyNA(y))
    stop("y should be a numeric vector without missing values.")
  if (missing(X) || !(inherits(X, "Matrix") || is.matrix(X)) || anyNA(X))
    stop("X should be a matrix without missing values.")
  if (missing(Sigma) || !(inherits(Sigma, "Matrix") || is.matrix(Sigma)) ||anyNA(Sigma))
    stop("Sigma should be a matrix without missing values.")
  if (!is.null(covs) && (!(inherits(covs, "Matrix") || is.matrix(covs)) || anyNA(covs)))
    stop("covs should be a numeric vector without missing values.")
  if (!is.numeric(nChunks) || length(nChunks) > 1 || nChunks != round(nChunks))
    stop("nChunks should be an integer")
  n <- length(y)
  ## Check dimensions.
  if (nrow(X) != n)
    stop("The number of elements in y should be identical to the number of rows in X")
  if (nrow(Sigma) != n || ncol(Sigma) != n)
    stop("The number of elements in y should be identical to the number of rows and columns in Sigma")
  if (!is.null(covs) && nrow(covs) != n)
    stop("The number of elements in y should be identical to the number of rows in covs")
  m <- ncol(X)
  ## If necessary convert input to Matrix
  if (is.matrix(X)) X <- as(X, "dgeMatrix")
  if (is.matrix(Sigma)) Sigma <- as(Sigma, "dsyMatrix")
  if (is.matrix(covs)) covs <- as(covs, "dgeMatrix")
  ## Number of chunks should be smaller than m.
  if (nChunks > m) nChunks <- ceiling(m / 2)
  fixCovs <- Matrix::cbind2(rep(1, n), covs)
  nCov <- ncol(fixCovs)
  M <- Matrix::solve(Matrix::chol(Sigma))
  ## Pre-multiply the phenotype (y) with t(M).
  tMy <- Matrix::crossprod(M, y)
  ## pre-multiply the intercept and covariates with t(M).
  tMfixCovs <- Matrix::crossprod(M, fixCovs)
  ## Pre-multiply the snp-matrix with t(M).
  tMX <- Matrix::crossprod(M, X)
  ## Matrix cookbook, 3.2.6 Rank-1 update of inverse of inner product.
  A <- Matrix::solve(Matrix::crossprod(tMfixCovs))
  vv <- Matrix::colSums(tMX ^ 2)
  vX <- Matrix::crossprod(tMfixCovs, tMX)
  nn <- 1 / (vv - Matrix::colSums(vX * (A %*% vX)))
  XtXinvLastRows <- Matrix::cbind2(- nn * Matrix::crossprod(vX, A), nn)
  Xty <- Matrix::cbind2(Matrix::Matrix(rep(as.numeric(Matrix::crossprod(tMfixCovs, tMy)),
    length(nn)), byrow = TRUE, ncol = nCov), Matrix::crossprod(tMX, tMy))
  betaVec <- Matrix::rowSums(XtXinvLastRows[, 1:nCov, drop = FALSE] * Xty[, 1:nCov, drop = FALSE]) +
    XtXinvLastRows[, 1 + nCov] * Xty[, 1 + nCov]
  ## Compute residuals and RSS over all markers.
  ResEnv <- lsfit(x = tMfixCovs, y = tMy, intercept = FALSE)$residuals
  RSSEnv <- sum(ResEnv ^ 2)
  ## QR decomposition of covariates.
  Q <- Matrix::qr.Q(Matrix::qr(tMfixCovs))
  tMQtQ <- Matrix::t(M %*% (Matrix::Diagonal(n = n) - Matrix::tcrossprod(Q)))
  ## Compute RSS per marker, breaking up X for speed.
  ## In case nChunks = 1 everything can be done in a single step. Otherwise loop over the chunks.
 if (nChunks > 1) {
    chunks <- split(1:m, c(rep(1:nChunks, each = m %/% nChunks),
      rep(nChunks, each = m %% nChunks)))
    ## In case nChunks = 1 everything can be done in a single step. Otherwise loop over the chunks.
    RSSFull <- unlist(lapply(seq_along(chunks), FUN = function(i) {
      tX <- tMQtQ %*% X[, chunks[[i]]]
      apply(tX, 2, function(x) {
        sum(lsfit(x = x, y = ResEnv, intercept = FALSE)$residuals ^ 2)})
    }))
  } else {
    tX <- tMQtQ %*% X
    RSSFull <- apply(tX, 2, function(x) {
      sum(lsfit(x = x, y = ResEnv, intercept = FALSE)$residuals ^ 2)})
  }
  ## Compute F and p values.
  df2 <- n - 1 - nCov
  FVal <- (RSSEnv - RSSFull) / RSSFull * df2
  pVal <- pf(q = FVal, df1 = 1, df2 = df2, lower.tail = FALSE)
  ## Compute R_LR^2 statistic from Sun et al 2010, heredity.
  RLR2  <- 1 - exp((RSSFull - RSSEnv) / n)
  ## Construct output data.frame.
  GLS <- data.frame(pValue = pVal,
    beta = betaVec,
    betaSe = sqrt(XtXinvLastRows[, 1 + nCov]),
    RLR2 = RLR2)
  return(GLS)
}

