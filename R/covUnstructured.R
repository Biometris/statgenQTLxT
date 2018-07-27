#' Compute unstructured covariance
#'
#' Compute unstructured covariance pairwise using \code{covPairwise} or using a
#' single model using \code{covUnstructured}.
#'
#' @inheritParams EMFA

#' @param X A covariate matrix, c being the number of covariates and n being the
#' number of genotypes.
#' @param fixDiag Should the diagonal of the covariate matrix be fixed during
#' calculations? -- NOT YET IMPLEMENTED
#' @param VeDiag Should Ve be a diagonale matrix?
#' @param corMat Should the output be a correlation matrix instead of a
#' covariance matrix?
#' @param parallel Should the computation of variance components be done in
#' parallel.
#'
#' @return A list of two matrices \code{Vg} and \code{Ve} containing genotypic
#' and environmental variance components respectively.
#'
#' @import utils
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @keywords internal
covUnstructured <- function(Y,
                            K,
                            X = NULL,
                            fixDiag = FALSE,
                            VeDiag = FALSE) {
  ## Check input.
  if (missing(Y) || !(is.matrix(Y) || inherits(Y, "Matrix"))) {
    stop("Y should be a matrix")
  }
  if (missing(K) || !(is.matrix(K) || inherits(K, "Matrix"))) {
    stop("K should be a matrix")
  }
  if (fixDiag) {
    warning("fixDiag = TRUE not implemented yet. Value set to FALSE")
    fixDiag <- FALSE
  }
  Y <- cbind(rownames(Y), as.data.frame(as.matrix(Y)), stringsAsFactors = FALSE)
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore Simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- cbind(rownames(X), as.data.frame(as.matrix(X)),
               stringsAsFactors = FALSE)
    dat <- merge(Y, X, by = "genotype")
  } else {
    dat <- Y
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(Y$genotype), unique(Y$genotype)]
  traits <- colnames(Y)[-1]
  nTrait <- length(traits)
  smpVar <- sapply(Y[-1], var)
  if (!is.null(X)) {
    ## Define formula for fixed part. ` needed to accommodate - in
    ## variable names.
    fixed <- as.formula(paste0("cbind(",
                               paste0(traits, collapse = ", "), ") ~ `",
                               paste(colnames(X)[-1], collapse = '` + `'), "`"))
  } else {
    fixed <- as.formula(paste0("cbind(", paste0(traits, collapse = ", "),
                               ") ~ 1"))
  }
  if (VeDiag) {
    rcov <- as.formula(~ diag(trait):units)
  } else {
    rcov <- as.formula(~ us(trait):units)
  }
  ## Fit model.
  sommerFit <- sommer::mmer2(fixed = fixed, random = ~ us(trait):g(genotype),
                             rcov = rcov, data = dat, G = list(genotype = K),
                             silent = TRUE, date.warning = FALSE)
  ## Extract components from fitted model.
  VgMat <- sommerFit$var.comp[[1]]
  VeMat <- sommerFit$var.comp[[2]]
  ## Keep diagonal for Vg and Ve away from 0.
  diag(VgMat)[diag(VgMat) <= 0] <- 1e-3 * smpVar[diag(VgMat) <= 0]
  diag(VeMat)[diag(VeMat) <= 0] <- 1e-3 * smpVar[diag(VeMat) <= 0]
  ## Make Vg and Ve positive definite
  VgMat <- Matrix::nearPD(VgMat)$mat
  VeMat <- Matrix::nearPD(VeMat)$mat
  colnames(VgMat) <- rownames(VgMat) <- traits
  colnames(VeMat) <- rownames(VeMat) <- traits
  return(list(Vg = VgMat, Ve = VeMat))
}

#' @rdname covUnstructured
#' @keywords internal
covPairwise <- function(Y,
                        K,
                        X = NULL,
                        fixDiag = FALSE,
                        corMat = FALSE,
                        parallel = FALSE) {
  ## Check input.
  if (missing(Y) || !(is.matrix(Y) || inherits(Y, "Matrix"))) {
    stop("Y should be a matrix")
  }
  if (missing(K) || !(is.matrix(K) || inherits(K, "Matrix"))) {
    stop("K should be a matrix")
  }
  if (fixDiag) {
    warning("fixDiag = TRUE not implemented yet. Value set to FALSE")
    fixDiag <- FALSE
  }
  `%op%` <- getOper(parallel && foreach::getDoParWorkers() > 1)
  Y <- cbind(rownames(Y), as.data.frame(as.matrix(Y)), stringsAsFactors = FALSE)
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore Simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- cbind(rownames(X), as.data.frame(as.matrix(X)),
               stringsAsFactors = FALSE)
    dat <- merge(Y, X, by = "genotype")
  } else {
    dat <- Y
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(Y$genotype), unique(Y$genotype)]
  traits <- colnames(Y)[-1]
  nTrait <- length(traits)
  smpVar <- sapply(Y[-1], var)
  VgVec <- VeVec <- vector(mode = "numeric", length = nTrait)
  for (i in 1:nTrait) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in varnames.
      fixed <- as.formula(paste(traits[i], " ~ `",
                                paste(colnames(X)[-1], collapse = '` + `'),
                                "`"))
    } else {
      fixed <- as.formula(paste(traits[i], " ~ 1"))
    }
    ## Fit model.
    sommerFit <- sommer::mmer2(fixed = fixed, random = ~ g(genotype),
                               data = dat, G = list(genotype = K),
                               silent = TRUE, date.warning = FALSE)
    ## Extract components from fitted model.
    VgVec[i] <- as.numeric(sommerFit$var.comp[[1]])
    VeVec[i] <- as.numeric(sommerFit$var.comp[[2]])
  }
  ## Keep diagonal for Vg and Ve away from 0.
  VgVec[VgVec <= 0] <- 1e-3 * smpVar[VgVec <= 0]
  VeVec[VeVec <= 0] <- 1e-3 * smpVar[VeVec <= 0]
  if (corMat) {
    ## Ones on the diagonal of resulting matrix.
    VgMat <- VeMat <- diag(x = 1, nrow = nTrait)
  } else {
    ## Computed values from univariate analysis on diagonal of resulting matrix.
    VgMat <- diag(x = VgVec)
    VeMat <- diag(x = VeVec)
  }
  rownames(VgMat) <- colnames(VgMat) <- traits
  rownames(VeMat) <- colnames(VeMat) <- traits
  ## For every combination of traits compute variance.
  modPW <- function(i, j) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in varnames.
      fixed <- formula(paste0("cbind(", traits[i], ", ", traits[j], ") ~ `",
                              paste(colnames(X)[-1], collapse = '` + `'), "`"))
    } else {
      fixed <- formula(paste0("cbind(", traits[i], ", ", traits[j], ") ~ 1"))
    }
    modFit <- sommer::mmer2(fixed = fixed, random = ~ us(trait):g(genotype),
                            rcov = ~ us(trait):units,
                            data = dat[, c("genotype", traits[i], traits[j],
                                           colnames(X)[-1])],
                            G = list(genotype = K),
                            init = list(diag(VgVec[c(i, j)]),
                                        diag(VeVec[c(i, j)])),
                            silent = TRUE, date.warning = FALSE)
    return(c(modFit$var.comp[[1]][1,2],
             modFit$var.comp[[2]][1,2]))
  }
  comb <- combn(x = seq_along(traits), m = 2)
  pwVar <- foreach::foreach(i = comb[1, ], j = comb[2, ],
                            .combine = "cbind", .packages = "sommer") %op% {
                              modPW(i, j)
                            }
  ## If there are only 2 traits pwVar is a vector, where it should be a matrix.
  if (nTrait == 2) {
    pwVar <- t(t(pwVar))
  }
  ## Fill VgMat using symmetry.
  VgMat[lower.tri(VgMat)] <- pwVar[1, ]
  VgMat[upper.tri(VgMat)] <- t(VgMat)[upper.tri(VgMat)]
  if (corMat) {
    VgMat <- cor(VgMat)
  }
  ## Fill VeMat using symmetry.
  VeMat[lower.tri(VeMat)] <- pwVar[2, ]
  VeMat[upper.tri(VeMat)] <- t(VeMat)[upper.tri(VeMat)]
  if (corMat) {
    VeMat <- cor(VeMat)
  }
  ## Make positive definite.
  VgMat <- Matrix::nearPD(VgMat, corr = corMat)$mat
  VeMat <- Matrix::nearPD(VeMat, corr = corMat)$mat
  return(list(Vg = VgMat, Ve = VeMat))
}







