#' Compute pairwise covariance
#'
#' Compute pairwise covariance.
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


## TO DO: univariate G-BLUPs; + correlations in case of non-convergence
## diagonal Ve
## p-values for correlations

#' @import utils

covPairwise <- function(Y,
  K,
  X = NULL,
  fixDiag = FALSE,
  corMat = FALSE,
  VeDiag = FALSE) {
  ## Check input.
  if (missing(Y) || !is.matrix(Y))
    stop("Y should be a matrix")
  if (missing(K) || !is.matrix(K))
    stop("K should be a matrix")
  if (fixDiag) {
    warning("fixDiag = TRUE not implemented yet. Value set to FALSE")
    fixDiag <- FALSE
  }
  Y <- tibble::rownames_to_column(as.data.frame(Y), var = "genotype")
  if (!is.null(X)) {
    X <- tibble::rownames_to_column(as.data.frame(X), var = "genotype")
    data <- merge(Y, X, by = "genotype")
  } else {
    data <- Y
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(Y$genotype), unique(Y$genotype)]
  traits <- colnames(Y)[-1]
  nTrait <- length(traits)
  convMat <- matrix(FALSE, nTrait, nTrait)
  colnames(convMat) <- rownames(convMat) <- traits
  if (!is.null(X)) {
    ## Define formula for fixed part. ` needed to accommodate - in variable names.
    fixed <- as.formula(paste0("cbind(", paste0(traits, collapse = ", "), ") ~ `",
      paste(colnames(X)[-1], collapse ='` + `'), "`"))
  } else {
    fixed <- as.formula(paste0("cbind(", paste0(traits, collapse = ", "), ") ~ 1"))
  }
  ## Fit model.
  sommerFit <- sommer::mmer2(fixed = fixed, random = ~ sommer::g(genotype),
    data = data, G = list(genotype = K), silent = TRUE)
  ## Extract components from fitted model.
  VgVec <- diag(sommerFit$var.comp[[1]])
  VeVec <- diag(sommerFit$var.comp[[2]])
  convVec <- sapply(X = traits, FUN = function(x) {sommerFit[[x]]$converge})
  if (corMat) {
    ## Ones on the diagonal of resulting matrix.
    VgMat <- VeMat <- diag(x = 1, nrow = nTrait)
  } else {
    ## Computed values from univariate analysis on diagonal of resulting matrix.
    VgMat <- diag(x = VgVec)
    VeMat <- diag(x = VeVec)
  }
  if (VeDiag) {
    rcov <- as.formula(~ diag(trait):units)
  } else {
    rcov <- as.formula(~ us(trait):units)
  }
  ## For every combination of traits compute variance.
  pwVar <- combn(
    traits, 2,
    FUN = function(idx) {
      if (!is.null(X)) {
        ## Define formula for fixed part. ` needed to accommodate - in variable names.
        fixed <- as.formula(paste0("cbind(", idx[[1]], ", ", idx[[2]], ") ~ `",
          paste(colnames(X)[-1], collapse ='` + `'), "`"))
      } else {
        fixed <- as.formula(paste0("cbind(", idx[[1]], ", ", idx[[2]], ") ~ 1"))
      }
      sommerFit <- sommer::mmer2(fixed = fixed, random = ~ us(trait):g(genotype),
        rcov = rcov, data = data,
        G = list(genotype = K), silent = TRUE)
      ## Extract components from fitted model.
      return(sommerFit$var.comp)
    }, simplify = FALSE)
  ## Fill VgMat using symmetry.
  VgMat[lower.tri(VgMat)] <- sapply(1:length(pwVar), FUN = function(x) {
    if (!corMat) pwVar[[x]][[1]][1, 2] else cov2cor(pwVar[[x]][[1]])[1, 2]
  })
  VgMat[upper.tri(VgMat)] <- t(VgMat)[upper.tri(VgMat)]
  ## Fill VeMat using symmetry.
  VeMat[lower.tri(VeMat)] <- sapply(1:length(pwVar), FUN = function(x) {
    if (!corMat) pwVar[[x]][[2]][1, 2] else cov2cor(pwVar[[x]][[2]])[1, 2]
  })
  VeMat[upper.tri(VeMat)] <- t(VeMat)[upper.tri(VeMat)]
  ## Make positive definite.
  VgMat <- as.matrix(Matrix::nearPD(VgMat, corr = corMat)$mat)
  VeMat <- as.matrix(Matrix::nearPD(VeMat, corr = corMat)$mat)
  ## Multiply by results from univariate analysis.
  VgMat <- tcrossprod(matrix(sqrt(VgVec))) * VgMat
  VeMat <- tcrossprod(matrix(sqrt(VeVec))) * VeMat
  ## Add row- and column names.
  colnames(VgMat) <- rownames(VgMat) <- traits
  colnames(VeMat) <- rownames(VeMat) <- traits
  return(list(Vg = VgMat, Ve = VeMat))
}



