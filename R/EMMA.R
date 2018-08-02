#' Compute REML estimates of variance components using EMMA algorithm.
#'
#' Using the EMMA algorithm as is Kang et al. (2008) compute REML estimates of
#' genetic and residual variance components.
#'
#' @param gData An object of class gData containing at least a data.frame
#' \code{pheno}. If \code{K} is not supplied a matrix \code{kinship} should be
#' in \code{gData}. If covariates are included then a data.frame covar is
#' needed as wel and if an extra snp is to be included as covariate (defined
#' in \code{snpName}) then a data.frame \code{markers} is also needed. Missing
#' values in \code{pheno} are allowed but will be excluded from the
#' calculations.
#' @param trait A trait for which to estimate variance components. This can be
#' either numeric index or character name of a column in \code{pheno}.
#' @param environment An environment for which to estimate variance components.
#' This can be either numeric index or character name of a list item in
#' \code{pheno}.
#' @param K An optional kinship matrix. If \code{NULL} then matrix
#' \code{kinship} in \code{gData} is used. If both \code{K} is provided and
#' \code{gData} contains a matrix \code{kinship} then \code{K} is used.
#' @param covar An optional vector of covariates taken into account when
#' estimating variance components. These can be either numeric indices or
#' character names of columns in \code{covar} in \code{gData}. If \code{NULL}
#' no covariates are used.
#' @param snpName An optional character string of a marker in \code{markers} in
#' \code{gData} to be included as covariate. If used the \code{gData} object
#' should contain a data.frame \code{markers}.
#' @param Z An optional incidence matrix mapping each observed phenotype to
#' one of inbred strains.
#' @param nGrids An integer indicating the number of intervals used for local
#' optimisation within the algorithm.
#' @param lLim A numerical value indicating the lower limit of the interval over
#' which optimisating is done.
#' @param uLim A numerical value indicating the upper limit of the interval over
#' which optimisating is done.
#' @param eps A numerical value used as computational tolerance in the
#' algorithm.
#'
#' @return A list with two components:
#' \itemize{
#' \item{\code{varcomp} a vector of genetic variance Vg and residual variance
#' Ve}
#' \item{\code{vcovMatrix} The variance covariance matrix corresponding to
#' the computed variances.}
#' }
#' @references Kang et al. (2008) (Efficient Control of Population Structure in
#' Model Organism Association Mapping. Genetics, March 2008, Vol. 178, no. 3,
#' p. 1709-1723
#' @import stats
#'
#' @keywords internal
EMMA <- function(gData,
                 trait,
                 environment,
                 K = NULL,
                 covar = NULL,
                 snpName = NULL,
                 Z = NULL,
                 nGrids = 100,
                 lLim = -10,
                 uLim = 10,
                 eps = 1e-10) {
  ## Checks.
  chkGData(gData, comps = "pheno")
  if (missing(environment) || length(environment) > 1 ||
      !(is.numeric(environment) || is.character(environment))) {
    stop("environment should be a single numeric or character.\n")
  }
  if ((is.character(environment) && !environment %in% names(gData$pheno)) ||
      (is.numeric(environment) && environment > length(gData$pheno))) {
    stop("environment should be a list item in pheno.\n")
  }
  if (missing(trait) || length(trait) > 1 || !(is.numeric(trait) ||
                                               is.character(trait))) {
    stop("trait should be a single numeric or character.\n")
  }
  if ((is.character(trait) && !trait %in%
       colnames(gData$pheno[[environment]])) ||
      (is.numeric(trait) && trait > ncol(gData$pheno[[environment]]))) {
    stop("trait should be a column in pheno.\n")
  }
  if (!is.null(K) && !(inherits(K, "Matrix") || is.matrix(K))) {
    stop("K should be a matrix.\n")
  }
  if (is.null(K) && is.null(gData$kinship)) {
    stop("gData contains no matrix kinship so K should be provided.\n")
  }
  if (!is.null(covar) && !is.numeric(covar) && !is.character(covar)) {
    stop("covar should be a numeric or character.\n")
  }
  if ((is.character(covar) && !all(covar %in% colnames(gData$covar))) ||
      (is.numeric(covar) && any(covar > ncol(gData$covar)))) {
    stop("covar should be columns in covar in gData.\n")
  }
  if (!is.null(snpName) && (length(snpName) > 1 || !is.character(snpName))) {
    stop("snpName should be a single character.\n")
  }
  if (!is.null(Z) && !is.matrix(Z)) {
    stop("Z should be a matrix.\n")
  }
  if (!is.null(nGrids) && (length(nGrids) > 1 || !is.numeric(nGrids) ||
                           nGrids != round(nGrids))) {
    stop("nGrids should be a single integer.\n")
  }
  if (!is.null(lLim) && (length(lLim) > 1 || !is.numeric(lLim))) {
    stop("lLim should be a single numerical value.\n")
  }
  if (!is.null(uLim) && (length(uLim) > 1 || !is.numeric(uLim))) {
    stop("uLim should be a single numerical value.\n")
  }
  if (lLim >= uLim) {
    stop("lLim should be smaller than uLim.\n")
  }
  if (!is.null(eps) && (length(eps) > 1 || !is.numeric(eps))) {
    stop("eps should be a single numerical value.\n")
  }
  ## Add column genotype to environment.
  phEnv <- gData$pheno[[environment]]
  ## Remove data with missings in trait or any of the covars.
  nonMiss <- phEnv$genotype[!is.na(phEnv[trait])]
  nonMissId <- which(!is.na(phEnv[trait]))
  if (!is.null(covar)) {
    misCov <- rownames(gData$covar)[rowSums(is.na(gData$covar[covar])) == 0]
    nonMiss <- nonMiss[nonMiss %in% misCov]
    nonMissId <- intersect(nonMissId, which(phEnv$genotype %in% misCov))
  }
  if (is.null(K)) {
    K <- gData$kinship[nonMiss, nonMiss]
  } else {
    K <- K[nonMiss, nonMiss]
  }
  y <- as(phEnv[nonMissId, trait], "dgeMatrix")
  ## Define intercept.
  X <- Matrix::Matrix(data = 1, nrow = length(nonMiss), ncol = 1)
  if (!is.null(covar)) {
    ## Add covars to intercept.
    X <- Matrix::cbind2(X, as(as.matrix(gData$covar[nonMiss, covar]),
                              "dgeMatrix"))
  }
  if (!is.null(snpName)) {
    ## Add extra snp to intercept + covars.
    X <- Matrix::cbind2(X, as.numeric(gData$markers[phEnv$genotype,
                                                    snpName][nonMiss]))
  }
  ## Check resulting X for singularity.
  if (!inherits(try(Matrix::solve(Matrix::crossprod(X)), silent = TRUE),
                "Matrix")) {
    warning("X is singular.")
    return(list(varcomp = c(0, 0), K = K))
  }
  ## Using notation of Kang et al.
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  ## Define intervals used for computing local optimums.
  m <- nGrids + 1
  logDelta <- (0:nGrids) / nGrids * (uLim - lLim) + lLim
  delta <- exp(logDelta)
  if (is.null(Z)) {
    ## Compute n-q non-zero eigenvalues and corresponding eigenvectors.
    eigR <- emmaEigenR(K = K, X = X)
    ## Define eta as in eqn. 7 of Kang.
    etas <- Matrix::crossprod(eigR$vectors, y)
    ## Define etas1 and etas2 to use same optimisating as for Z not NULL.
    etas1 <- etas
    etas2 <- 0
    ## Compute square of eta for usage in vectorised form.
    etasQ <- etas ^ 2
    ## Define lambda as in eqn. 7 of Kang for usage in vectorised form.
    lambdas <- Matrix::Matrix(eigR$values, nrow = n - q, ncol = m) +
      Matrix::Matrix(delta, nrow = n - q, ncol = m, byrow = TRUE)
    ## Compute derivative of LL as in eqn. 9 of Kang for all grid endpoints.
    dLL <- 0.5 * delta * ((n - q) * Matrix::colSums(etasQ / lambdas ^ 2) /
                            Matrix::colSums(etasQ / lambdas) -
                            Matrix::colSums(1 / lambdas))
  } else {
    ## Compute n-q non-zero eigenvalues and corresponding eigenvectors.
    eigR <- emmaEigenRZ(Z = Z, K = K, X = X)
    ## Define eta as in eqn. 7 of Kang.
    etas <- Matrix::crossprod(eigR$vectors, y)
    ## Split etas
    etas1 <- etas[1:(t - q)]
    etas2 <- sum(etas[(t - q + 1):(n - q)] ^ 2)
    ## Compute square of eta for usage in vectorised form.
    etasQ <- Matrix::Matrix(etas1 ^ 2, nrow = t - q, ncol = m)
    ## Define lambda as in eqn. 7 of Kang for usage in vectorised form.
    lambdas <- Matrix::Matrix(eigR$values, nrow = t - q, ncol = m) +
      Matrix::Matrix(delta, nrow = t - q, ncol = m, byrow = TRUE)
    ## Compute derivative of LL as in eqn. 9 of Kang for all grid endpoints.
    dLL <- 0.5 * delta * ((n - q) *
                            (Matrix::colSums(etasQ / (lambdas ^ 2)) + etas2 /
                               (delta ^ 2)) /
                            (Matrix::colSums(etasQ / lambdas) + etas2 / delta) -
                            (Matrix::colSums(1 / lambdas) + (n - t) / delta))
  }
  ## Find optimum of LL
  optLogDelta <- numeric(0)
  optLL <- numeric(0)
  ## Check first item in dLL. If < eps include LL value as possible optimum.
  if (dLL[1] < eps) {
    optLogDelta <- c(optLogDelta, lLim)
    optLL <- c(optLL, emmaREMLLL(logDelta = lLim, lambda = eigR$values,
                                 etas1 = etas1, n = n, t = t, etas2 = etas2))
  }
  ## Check last item in dLL. If > - eps include LL value as possible optimum.
  if (dLL[m] > -eps) {
    optLogDelta <- c(optLogDelta, uLim)
    optLL <- c(optLL, emmaREMLLL(logDelta = uLim, lambda = eigR$values,
                                 etas1 = etas, n = 0, t = 0, etas2 = 0))
  }
  ## If derivative changes sign on an interval, compute local optimum for LL
  ## on that interval and add it to possible optima.
  for (i in 1:(m - 1)) {
    if ((dLL[i] > 0 && dLL[i + 1] < 0) || dLL[i] * dLL[i + 1] < eps ^ 2) {
      r <- optimise(f = emmaREMLLL, lower = logDelta[i],
                    upper = logDelta[i + 1], lambda = eigR$values,
                    etas1 = etas, n = 0, t = 0, etas2 = 0, maximum = TRUE)
      optLogDelta <- c(optLogDelta, r$maximum)
      optLL <- c(optLL, r$objective)
    }
  }
  ## Compute absolute LL maximum from possible optima.
  maxDelta <- exp(optLogDelta[which.max(optLL)])
  maxLL <- max(optLL)
  ## Compute variance components.
  if (is.null(Z)) {
    maxVg <- sum(etas ^ 2 / (eigR$values + maxDelta)) / (n - q)
  } else {
    maxVg <- (sum(etas1 ^ 2 / (eigR$values + maxDelta)) + etas2 / maxDelta) /
      (n - q)
  }
  maxVe <- maxVg * maxDelta
  vcovMatrix <- maxVg * K + Matrix::Diagonal(n = nrow(K), x = maxVe)
  rownames(vcovMatrix) <- colnames(vcovMatrix) <- rownames(K)
  return(list(varComp = c(Vg = maxVg, Ve = maxVe), vcovMatrix = vcovMatrix))
}

#' EMMA helper functions
#'
#' Helper functions for computing REML estimates of genetic and residual
#' variance components using the EMMA algorithm.
#'
#' @inheritParams EMMA
#' @param X a q x n covariate matrix, q being the number of covariates and n
#' being the number of genotypes. q has to be at least one (typically an
#' intercept).
#'
#' @keywords internal
emmaEigenR <- function(K,
                       X) {
  n <- nrow(X)
  q <- ncol(X)
  ## Compute n-q non-zero eigenvalues of SHS as defined in eqn. 5 of Kang.
  if (q == 1) {
    S <- Matrix::Diagonal(n = n) - 1 / n
  } else {
    S <- Matrix::Diagonal(n = n) - X %*%
      Matrix::solve(Matrix::crossprod(X), Matrix::t(X))
  }
  eig <- eigen(S %*% (K + Matrix::Diagonal(n = n)) %*% S, symmetric = TRUE)
  if (is.complex(eig$values))
    stop("Complex eigen values found.\n")
  return(list(values = eig$values[1:(n - q)] - 1,
              vectors = eig$vectors[, 1:(n - q)]))
}

#' @keywords internal
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
              vectors = qr.Q(qr(cbind(SZ %*% eig$vectors[, 1:(t - q)], qrX)),
                             complete = TRUE)[, c(1:(t - q), (t + 1):n)]))
}

emmaREMLLL <- function(logDelta, lambda, etas1, n, t, etas2) {
  ## Compute the REML LL as in eqn. 7 of Kang.
  nq <- length(etas1) + n - t
  delta <- exp(logDelta)
  lDelta <- lambda + delta
  return(0.5 * (nq * (log(nq / (2 * pi)) - 1 -
                        log(sum(etas1 ^ 2 / (lDelta)) + etas2 / delta)) -
                  (sum(log(lDelta)) + (n - t) * logDelta)))
}

