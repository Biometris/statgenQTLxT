#' Compute unstructured covariance
#'
#' Compute unstructured covariance pairwise using \code{covPairwise} or using a
#' single model using \code{covUnstr}.
#'
#' @param Y An n x p matrix of observed phenotypes, on p traits or trials
#' for n individuals. No missing values are allowed.
#' @param K An n x n kinship matrix.
#' @param X An n x c covariate matrix, c being the number of covariates and n
#' being the number of genotypes.
#' @param fixDiag Should the diagonal of the covariate matrix be fixed during
#' calculations? -- NOT YET IMPLEMENTED
#' @param VeDiag Should Ve be a diagonal matrix?
#' @param corMat Should the output be a correlation matrix instead of a
#' covariance matrix?
#' @param parallel Should the computation of variance components be done in
#' parallel?
#'
#' @return A list of two matrices \code{Vg} and \code{Ve} containing genotypic
#' and environmental variance components respectively.
#'
#' @import utils stats
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @keywords internal
covUnstr <- function(Y,
                     K,
                     X = NULL,
                     fixDiag = FALSE,
                     VeDiag = FALSE) {
  ## Check input.
  if (missing(Y) || !(is.matrix(Y) || inherits(Y, "matrix"))) {
    stop("Y should be a matrix")
  }
  if (missing(K) || !(is.matrix(K) || inherits(K, "matrix"))) {
    stop("K should be a matrix")
  }
  if (fixDiag) {
    warning("fixDiag = TRUE not implemented yet. Value set to FALSE")
    fixDiag <- FALSE
  }
  ## Add genotype to data to be used in fitting model.
  ## This silently converts Y to a data.frame.
  dat <- cbind(genotype = rownames(Y), as.data.frame(as.matrix(Y)),
               stringsAsFactors = FALSE)
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore Simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- cbind(genotype = rownames(X), as.data.frame(as.matrix(X)),
               stringsAsFactors = FALSE)
    dat <- merge(dat, X, by = "genotype")
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(dat$genotype), unique(dat$genotype)]
  traits <- colnames(Y)
  nTrait <- length(traits)
  smpVar <- apply(X = Y, MARGIN = 2, FUN = var)
  ratio <- vector(mode = "numeric", length = nTrait)
  for (i in 1:nTrait) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in varnames.
      fixed <- formula(paste0(traits[i], " ~ `",
                              paste(colnames(X)[-1], collapse = '` + `'),
                              "`"))
    } else {
      fixed <- formula(paste(traits[i], " ~ 1"))
    }
    ## Fit model.
    modFit <- sommer::mmer(fixed = fixed,
                           random = ~sommer::vsr(genotype, Gu = as.matrix(K)),
                           data = dat, verbose = FALSE, date.warning = FALSE)
    ## Extract components from fitted model.
    Vg <- as.numeric(modFit$sigma[[1]])
    Ve <- as.numeric(modFit$sigma[[2]])
    ratio[i] <- Vg / (Vg + Ve)
  }
  if (all(ratio < 0.05)) {
    Vg <- matrix(data = 0, nrow = nTrait, ncol = nTrait)
    Ve <- cov(Y)
  } else {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in
      ## variable names.
      fixed <- formula(paste0("cbind(",
                              paste0(traits, collapse = ", "), ") ~ `",
                              paste(colnames(X)[-1], collapse = '` + `'), "`"))
    } else {
      fixed <- formula(paste0("cbind(", paste0(traits, collapse = ", "),
                              ") ~ 1"))
    }
    if (VeDiag) {
      rcov <- formula(paste0("~sommer::vsr(units, Gtc = diag(", nTrait, "))"))
    } else {
      rcov <- formula(paste0("~sommer::vsr(units,
                              Gtc = sommer::unsm(", nTrait, "))"))
    }
    random <- formula(paste0("~sommer::vsr(genotype, Gu = as.matrix(K),
                              Gtc = sommer::unsm(", nTrait, "))"))
    ## Fit model.
    modFit <- sommer::mmer(fixed = fixed, random = random, rcov = rcov,
                           data = dat, verbose = FALSE, dateWarning = FALSE,
                           ## This is not really a good idea, but for now it
                           ## is better than nothing.
                           ## Has to be solved in sommer.
                           tolParInv = 1e-6, method = "AI")
    ## Extract components from fitted model.
    VgMat <- modFit$sigma[[1]]
    VeMat <- modFit$sigma[[2]]
  }
  ## Keep diagonal for Vg and Ve away from 0.
  diag(VgMat)[diag(VgMat) <= 0] <- 1e-6 * smpVar[diag(VgMat) <= 0]
  diag(VeMat)[diag(VeMat) <= 0] <- 1e-6 * smpVar[diag(VeMat) <= 0]
  ## Make Vg and Ve positive definite.
  VgMat <- nearestPD(VgMat)
  VeMat <- nearestPD(VeMat)
  colnames(VgMat) <- rownames(VgMat) <- traits
  colnames(VeMat) <- rownames(VeMat) <- traits
  return(list(Vg = VgMat, Ve = VeMat))
}

#' @rdname covUnstr
#'
#' @importFrom stats formula
#' @keywords internal
covPW <- function(Y,
                  K,
                  X = NULL,
                  fixDiag = FALSE,
                  corMat = FALSE,
                  parallel = FALSE) {
  ## Check input.
  if (missing(Y) || !(is.matrix(Y) || inherits(Y, "matrix"))) {
    stop("Y should be a matrix")
  }
  if (missing(K) || !(is.matrix(K) || inherits(K, "matrix"))) {
    stop("K should be a matrix")
  }
  if (fixDiag) {
    warning("fixDiag = TRUE not implemented yet. Value set to FALSE")
    fixDiag <- FALSE
  }
  `%op%` <- getOper(parallel && foreach::getDoParWorkers() > 1)
  dat <- cbind(genotype = rownames(Y), as.data.frame(as.matrix(Y)),
               stringsAsFactors = FALSE)
  if (!is.null(X)) {
    ## sommer cannot handle column names with special characters.
    ## Therefore simplify column names in X.
    colnames(X) <- make.names(colnames(X), unique = TRUE)
    X <- cbind(genotype = rownames(X), as.data.frame(as.matrix(X)),
               stringsAsFactors = FALSE)
    dat <- merge(dat, X, by = "genotype")
  }
  ## Restrict K to genotypes in Y.
  K <- K[unique(dat$genotype), unique(dat$genotype)]
  traits <- colnames(Y)
  nTrait <- length(traits)
  smpVar <- apply(X = Y, MARGIN = 2, FUN = var)
  VgVec <- VeVec <- vector(mode = "numeric", length = nTrait)
  for (i in 1:nTrait) {
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in varnames.
      fixed <- formula(paste0(traits[i], " ~ `",
                              paste(colnames(X)[-1], collapse = '` + `'),
                              "`"))
    } else {
      fixed <- formula(paste(traits[i], " ~ 1"))
    }
    ## Fit model.
    modFit <- sommer::mmer(fixed = fixed,
                           random = ~sommer::vsr(genotype, Gu = as.matrix(K)),
                           data = dat, verbose = FALSE, date.warning = FALSE)
    ## Extract components from fitted model.
    VgVec[i] <- as.numeric(modFit$sigma[[1]])
    VeVec[i] <- as.numeric(modFit$sigma[[2]])
  }
  ## Keep diagonal for Vg and Ve away from 0.
  VgVec[VgVec <= 0] <- 1e-6 * smpVar[VgVec <= 0]
  VeVec[VeVec <= 0] <- 1e-6 * smpVar[VeVec <= 0]
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
    tolParInv <- 1e-6
    cat(i, " ", j, "\n")
    if (!is.null(X)) {
      ## Define formula for fixed part. ` needed to accommodate - in varnames.
      fixed <- formula(paste0("cbind(", traits[i], ", ", traits[j], ") ~ `",
                              paste(colnames(X)[-1], collapse = '` + `'), "`"))
    } else {
      fixed <- formula(paste0("cbind(", traits[i], ", ", traits[j], ") ~ 1"))
    }
    while (TRUE) {
      ## Sometimes sommer models won't fit with the default value for
      ## tolParInv. It is then suggested to increase the value of this
      ## parameter to assure a fit. We don't want a fixed high value of
      ## tolParInv, so we increase it in steps if needed by taking the sqrt.
      msg <- capture.output(
        modFit <-
          sommer::mmer(fixed = fixed,
                       random = ~ sommer::vsr(genotype, Gu = as.matrix(K),
                                              Gtc = sommer::unsm(2),
                                              Gti = VgMat[c(i, j), c(i, j)]),
                       rcov = ~ sommer::vsr(units, Gtc = sommer::unsm(2),
                                            Gti = VeMat[c(i, j), c(i, j)]),
                       data = dat[, c("genotype", traits[i], traits[j],
                                      colnames(X)[-1])],
                       verbose = FALSE, dateWarning = FALSE,
                       tolParInv = tolParInv)
      )
      if (length(modFit) == 0) {
        tolParInv <- sqrt(tolParInv)
      } else {
        break
      }
    }
    return(c(modFit$sigma[[1]][1, 2],
             modFit$sigma[[2]][1, 2]))
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
  VgMat <- nearestPD(VgMat, corr = corMat)
  VeMat <- nearestPD(VeMat, corr = corMat)
  colnames(VgMat) <- rownames(VgMat) <- traits
  colnames(VeMat) <- rownames(VeMat) <- traits
  return(list(Vg = VgMat, Ve = VeMat))
}

