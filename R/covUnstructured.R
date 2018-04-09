#' Compute unstructured covariance
#'
#' Compute unstructured covariance pairwise using \code{covPairwise} or using a single model using
#' \code{covUnstructured}.
#'
#' @inheritParams EMFA

#' @param X a covariate matrix, c being the number of covariates and n being the number
#' of genotypes.
#' @param fixDiag should the diagonal of the covariate matrix be fixed during calculations?
#' -- NOT YET IMPLEMENTED
#' @param VeDiag should Ve be a diagonale matrix?
#' @param corMat should the output be a correlation matrix instead of a covariance matrix?
#'
#' @return a list of two matrices \code{Vg} and \code{Ve} containing genotypic and environmental
#' variance components respectively.
#'
#' @import utils
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
  Y <- tibble::rownames_to_column(as.data.frame(as.matrix(Y)), var = "genotype")
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore Simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- tibble::rownames_to_column(as.data.frame(as.matrix(X)), var = "genotype")
    data <- merge(Y, X, by = "genotype")
  } else {
    data <- Y
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(Y$genotype), unique(Y$genotype)]
  traits <- colnames(Y)[-1]
  nTrait <- length(traits)
  smpVar <- sapply(Y[-1], var)
  if (!is.null(X)) {
    ## Define formula for fixed part. ` needed to accommodate - in variable names.
    fixed <- as.formula(paste0("cbind(", paste0(traits, collapse = ", "), ") ~ `",
                               paste(colnames(X)[-1], collapse = '` + `'), "`"))
  } else {
    fixed <- as.formula(paste0("cbind(", paste0(traits, collapse = ", "), ") ~ 1"))
  }
  if (VeDiag) {
    rcov <- as.formula(~ diag(trait):units)
  } else {
    rcov <- as.formula(~ us(trait):units)
  }
  ## Fit model.
  sommerFit <- sommer::mmer2(fixed = fixed, random = ~ us(trait):g(genotype),
                             rcov = rcov, data = data, G = list(genotype = K),
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
                        corMat = FALSE) {
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
  Y <- tibble::rownames_to_column(as.data.frame(as.matrix(Y)), var = "genotype")
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore Simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- tibble::rownames_to_column(as.data.frame(as.matrix(X)), var = "genotype")
    data <- merge(Y, X, by = "genotype")
  } else {
    data <- Y
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(Y$genotype), unique(Y$genotype)]
  traits <- colnames(Y)[-1]
  nTrait <- length(traits)
  smpVar <- sapply(Y[-1], var)
  VgVec <- VeVec <- vector(mode = "numeric", length = nTrait)
  for (i in 1:nTrait) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in variable names.
      fixed <- as.formula(paste(traits[i], " ~ `",
                                paste(colnames(X)[-1], collapse = '` + `'),
                                "`"))
    } else {
      fixed <- as.formula(paste(traits[i], " ~ 1"))
    }
    ## Fit model.
    sommerFit <- sommer::mmer2(fixed = fixed, random = ~ g(genotype),
                               data = data, G = list(genotype = K),
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
  pwVar <- combn(x = traits, m = 2, FUN = function(idx) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in variable names.
      fixed <- as.formula(paste0("cbind(", idx[[1]], ", ", idx[[2]], ") ~ `",
                                 paste(colnames(X)[-1], collapse = '` + `'), "`"))
    } else {
      fixed <- as.formula(paste0("cbind(", idx[[1]], ", ", idx[[2]], ") ~ 1"))
    }
    sommerFit <- sommer::mmer2(fixed = fixed, random = ~ us(trait):g(genotype),
                               rcov = ~ us(trait):units, data = data,
                               G = list(genotype = K), silent = TRUE,
                               init = list(diag(diag(VgMat[idx, idx])),
                                           diag(diag(VeMat[idx, idx]))),
                               date.warning = FALSE)
    ## Extract components from fitted model.
    return(sommerFit$var.comp)
  }, simplify = FALSE)
  ## Fill VgMat using symmetry.
  VgMat[lower.tri(VgMat)] <- sapply(1:length(pwVar), FUN = function(x) {
    if (!corMat) pwVar[[x]][[1]][1, 2] else {
      if (pwVar[[x]][[1]][1, 2] == 0) {
        0
      } else {
        Matrix::cov2cor(Matrix::nearPD(pwVar[[x]][[1]])$mat)[1, 2]
      }
    }
  })
  VgMat[upper.tri(VgMat)] <- t(VgMat)[upper.tri(VgMat)]
  ## Fill VeMat using symmetry.
  VeMat[lower.tri(VeMat)] <- sapply(1:length(pwVar), FUN = function(x) {
    if (!corMat) pwVar[[x]][[2]][1, 2] else {
      if (pwVar[[x]][[2]][1, 2] == 0) {
        0
      } else {
        Matrix::cov2cor(Matrix::nearPD(pwVar[[x]][[2]])$mat)[1, 2]
      }
    }
  })
  VeMat[upper.tri(VeMat)] <- t(VeMat)[upper.tri(VeMat)]
  ## Make positive definite.
  VgMat <- Matrix::nearPD(VgMat, corr = corMat)$mat
  VeMat <- Matrix::nearPD(VeMat, corr = corMat)$mat
  return(list(Vg = VgMat, Ve = VeMat))
}



