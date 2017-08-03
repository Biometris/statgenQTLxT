#' Compute REML estimates of variance components using EMMA algorithm.
#'
#' Using the EMMA algorithm as is Kang et al. compute REML estimates of genetic and residual
#' variance components.
#'
#' @param gData an object of class gData containing at least a data.frame \code{pheno}. If \code{K} is
#' not supplied a matrix \code{kinship} should be in \code{gData}. If covariates are included
#' then a data.frame covar is needed as wel and if an extra snp is to be included as covariate
#' (defined in \code{snpName}) then a data.frame \code{markers} is also needed. Missing values
#' in \code{pheno} are allowed but will be excluded from the calculations.
#' @param trait a trait for which to estimate variance components. This can be either numeric index
#' or character name of a column in \code{pheno}.
#' @param field a field for which to estimate variance components. This can be either numeric index
#' or character name of a list item in \code{pheno}.
#' @param K an optional kinship matrix. If \code{NULL} then matrix \code{kinship} in \code{gData} is used.
#' If both \code{K} is provided and \code{gData} contains a matrix \code{kinship} then \code{K} is used.
#' @param covar an optional vector of covariates taken into account when estimating variance components.
#' These can be either numeric indices or character names of columns in \code{covar} in \code{gData}.
#' If \code{NULL} no covariates are used.
#' @param snpName an optional character name of a marker in \code{markers} in \code{gData} to be
#' included as covariate. If used the \code{gData} object should contain a data.frame \code{markers}.
#' @param Z an optional incidence matrix mapping each observed phenotype to one of inbred strains.
#' @param nGrids an integer indicating the number of intervals used for local optimisation within the
#' algorithm.
#' @param lLim a numeric value indicating the lower limit of the interval over which optimisating is
#' done.
#' @param uLim a numeric value indicating the upper limit of the interval over which optimisating is
#' done.
#' @param eps a numeric value used as computational tolerance in the algorithm.
#'
#' @return a list with two components:
#' \itemize{
#' \item{\code{varcomp} a vector of genetic variance Vg and residual variance Ve}
#' \item{\code{K} the kinship matrix used in the algorithm (entries for missing values in the
#' original kinship matrix are filtered out.)}
#' }
#' @references Kang et al. (2008) (Efficient Control of Population Structure in Model Organism
#' Association Mapping. Genetics, March 2008, Vol. 178, no. 3, p. 1709-1723

#' @import stats

# trait <- "anthesis.SYNHUN_2013_drought"
# snpName <- "AX-90549756"

runEmma <- function(gData,
  trait,
  field,
  K = NULL,
  covar = NULL,
  snpName = NULL,
  Z = NULL,
  nGrids = 100,
  lLim = - 10,
  uLim = 10,
  eps = 1e-10) {

  ## Check input
  if(missing(gData) || !is.gData(gData) || is.null(gData$pheno))
    stop("gData should be a valid gData object with at least pheno included.\n")
  if(missing(field) || length(field) > 1 || !(is.numeric(field) || is.character(field)))
    stop("field should be a single numeric or character.\n")
  if ((is.character(field) && !field %in% names(gData$pheno)) ||
      (is.numeric(field) && field > length(gData$pheno)))
    stop("field should be a list item in pheno.\n")
  if(missing(trait) || length(trait) > 1 || !(is.numeric(trait) || is.character(trait)))
    stop("trait should be a single numeric or character.\n")
  if ((is.character(trait) && !trait %in% colnames(gData$pheno[[field]])) ||
      (is.numeric(trait) && trait > ncol(gData$pheno[[field]])))
    stop("trait should be a column in pheno.\n")
  if (!is.null(K) && !is.matrix(K))
    stop("K should be a matrix.\n")
  if (is.null(K) && is.null(gData$kinship))
    stop("gData contains no matrix kinship so K should be provided.\n")
  if(!is.null(covar) && !is.numeric(covar) && !is.character(covar))
    stop("covar should be a numeric or character.\n")
  if ((is.character(covar) && !all(covar %in% colnames(covar))) ||
      (is.numeric(covar) && any(covar > ncol(covar))))
    stop("covar should be a columns in covar in gData.\n")
  if(!is.null(snpName) && (length(snpName) > 1 || !is.character(snpName)))
    stop("snpName should be a single character.\n")
  if(!is.null(Z) && !is.matrix(Z))
    stop("Z should be a matrix.\n")
  if(!is.null(nGrids) && (length(nGrids) > 1 || !is.numeric(nGrids) || nGrids != round(nGrids)))
    stop("nGrids should be a single integer.\n")
  if(!is.null(lLim) && (length(lLim) > 1 || !is.numeric(lLim)))
    stop("lLim should be a single numeric value.\n")
  if(!is.null(uLim) && (length(uLim) > 1 || !is.numeric(uLim)))
    stop("uLim should be a single numeric value.\n")
  if (lLim >= uLim)
    stop("lLim should be smaller than uLim.\n")
  if(!is.null(eps) && (length(eps) > 1 || !is.numeric(eps)))
    stop("eps should be a single numeric value.\n")

  ## Add column genotype to field.
  phenoField <- tibble::rownames_to_column(as.data.frame(gData$pheno[[field]]), var = "genotype")

  ## Remove data with missings in trait or any of the covars.
  nonMissing <- which(!is.na(phenoField[trait]))
  if (!is.null(covar)) {
    nonMissing <- intersect(nonMissing, which(rowSums(is.na(gData$covar[covar])) == 0))
  }
  if (is.null(K)) {
    K <- gData$kinship[nonMissing, nonMissing]
  } else {
    K <- K[nonMissing, nonMissing]
  }

  y <- phenoField[nonMissing, trait]

  ## Define intercept.
  X <- rep(1, length(nonMissing))
  if (!is.null(covar)) {
    ## Add covars to intercept.
    X <- cbind(X, gData$covar[nonMissing, covar])
  }
  if (!is.null(snpName)) {
    ## Add extra snp to intercept + covars.
    X <- cbind(X, as.numeric(gData$markers[phenoField$genotype, snpName][nonMissing]))
  }
  X <- as.matrix(X)
  ## Check resulting X for singularity.
  if (class(try(solve(crossprod(X)), silent = TRUE)) != "matrix") {
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
    etas <- crossprod(eigR$vectors, y)
    ## Define etas1 and etas2 to use same optimisating as for Z not NULL.
    etas1 <- etas
    etas2 <- 0
    ## Compute square of eta for usage in vectorised form.
    etasQ <- matrix(etas ^ 2 , nrow = n - q, ncol = m)
    ## Define lambda as in eqn. 7 of Kang for usage in vectorised form.
    lambdas <- matrix(eigR$values, nrow = n - q, ncol = m) +
      matrix(delta, nrow = n - q, ncol = m, byrow = TRUE)
    ## Compute derivative of LL as in eqn. 9 of Kang for all grid endpoints.
    dLL <- 0.5 * delta * ((n - q) * colSums(etasQ / lambdas ^ 2) / colSums(etasQ / lambdas) - colSums(1 / lambdas))
  } else {
    ## Compute n-q non-zero eigenvalues and corresponding eigenvectors.
    eigR <- emmaEigenRZ(Z = Z, K = K, X = X)
    ## Define eta as in eqn. 7 of Kang.
    etas <- crossprod(eigR$vectors, y)
    ## Split etas
    etas1 <- etas[1:(t - q)]
    etas2 <- sum(etas[(t - q + 1):(n - q)] ^ 2)
    ## Compute square of eta for usage in vectorised form.
    etasQ <- matrix(etas1 ^ 2, nrow = t - q, ncol = m)
    ## Define lambda as in eqn. 7 of Kang for usage in vectorised form.
    lambdas <- matrix(eigR$values, nrow = t - q, ncol = m) +
      matrix(delta, nrow = t - q, ncol = m, byrow = TRUE)
    ## Compute derivative of LL as in eqn. 9 of Kang for all grid endpoints.
    dLL <- 0.5 * delta * ((n - q) * (colSums(etasQ / (lambdas ^ 2)) + etas2 / (delta ^ 2)) /
        (colSums(etasQ / lambdas) + etas2 / delta) - (colSums(1 / lambdas) + (n - t) / delta))
  }
  ## Find optimum of LL
  optLogDelta <- numeric(0)
  optLL <- numeric(0)
  ## Check first item in dLL. If < eps include LL value as possible optima.
  if (dLL[1] < eps) {
    optLogDelta <- c(optLogDelta, lLim)
    optLL <- c(optLL, emmaREMLLL(logDelta = lLim, lambda = eigR$values, etas1 = etas1,
      n = n, t = t, etas2 = etas2))
  }
  ## Check last item in dLL. If > - eps include LL value as possible optima.
  if (dLL[m] > - eps) {
    optLogDelta <- c(optLogDelta, uLim)
    optLL <- c(optLL, emmaREMLLL(logDelta = uLim, lambda = eigR$values, etas1 = etas,
      n = 0, t = 0, etas2 = 0))
  }
  ## If derivative changes sign on an interval, compute local optimum for LL on that
  ## interval and add it to possible optima.
  for(i in 1:(m - 1)) {
    if (dLL[i] > 0 && dLL[i + 1] < 0 && dLL[i] * dLL[i + 1] < - eps ^ 2) {
      r <- optimise(emmaREMLLL, lower = logDelta[i], upper = logDelta[i + 1],
        lambda = eigR$values, etas1 = etas, n = 0, t = 0, etas2 = 0, maximum = TRUE)
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
    maxVg <- (sum(etas1 ^ 2 /(eigR$values + maxDelta)) + etas2 / maxDelta) / (n - q)
  }
  maxVe <- maxVg * maxDelta

  return(list(varcomp = c(Vg = maxVg, Ve = maxVe), K = K))
}














