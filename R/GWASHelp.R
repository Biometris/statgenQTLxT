#' Estimate variance components in multi trait GWAS
#'
#' Helper function for estimating variance components in multi trait GWAS.
#'
#' @noRd
#' @keywords internal
estVarComp <- function(GLSMethod,
                       covModel,
                       Y,
                       K,
                       X,
                       VeDiag,
                       snpCov,
                       XRed,
                       parallel,
                       maxIter,
                       mG,
                       mE,
                       chrs) {
  if (is.null(snpCov)) {
    varCompRed <- NULL
  }
  if (covModel %in% c("unst", "pw")) {
    ## Sommer always adds an intercept so remove it from X.
    XSom <- if (ncol(X) == 1) NULL else X[, -1, drop = FALSE]
    if (!is.null(snpCov)) {
      XRedSom <- if (ncol(XRed) == 1) NULL else XRed[, -1, drop = FALSE]
    }
  }
  if (covModel == "unst") {
    traitsMod <- restrictTraits(Y = Y, X = XSom, K = K)
    Y <- Y[,  traitsMod]
    ## Unstructured models.
    if (GLSMethod == "single") {
      varComp <- covUnstr(Y = Y, K = K, X = XSom, fixDiag = FALSE,
                          VeDiag = VeDiag)
      if (!is.null(snpCov)) {
        varCompRed <- covUnstr(Y = Y, K = K, X = XRedSom, fixDiag = FALSE,
                               VeDiag = VeDiag)
      }
    } else if (GLSMethod == "multi") {
      ## Unstructured models.
      varComp <- sapply(X = chrs, FUN = function(chr) {
        covUnstr(Y = Y, K = K[[which(chrs == chr)]],
                 X = XSom, fixDiag = FALSE, VeDiag = VeDiag)
      }, simplify = FALSE)
      if (!is.null(snpCov)) {
        varCompRed <- sapply(X = chrs, FUN = function(chr) {
          covUnstr(Y = Y, K = K[[which(chrs == chr)]],
                   X = XRedSom, fixDiag = FALSE, VeDiag = VeDiag)
        }, simplify = FALSE)
      }
    }
  } else if (covModel == "pw") {
    ## Unstructured (pairwise) models.
    if (GLSMethod == "single") {
      varComp <- covPW(Y = Y, K = K, X = XSom, fixDiag = FALSE, corMat = FALSE,
                       parallel = parallel)
      if (!is.null(snpCov)) {
        varCompRed <- covPW(Y = Y, K = K, X = XRedSom, fixDiag = FALSE,
                            corMat = FALSE, parallel = parallel)
      }
    } else if (GLSMethod == "multi") {
      varComp <- sapply(X = chrs, FUN = function(chr) {
        covPW(Y = Y, K = K[[which(chrs == chr)]], X = XSom, fixDiag = FALSE,
              corMat = FALSE, parallel = parallel)
      }, simplify = FALSE)
      if (!is.null(snpCov)) {
        varCompRed <- sapply(X = chrs, FUN = function(chr) {
          covPW(Y = Y, K = K[[which(chrs == chr)]], XRedSom, fixDiag = FALSE,
                corMat = FALSE, parallel = parallel)
        }, simplify = FALSE)
      }
    }
  } else if (covModel == "fa") {
    ## FA models.
    maxDiag <- 1e6 * max(abs(solve(var(Y))))
    if (GLSMethod == "single") {
      ## Including snpCovariates.
      varComp <- EMFA(y = Y, k = K, size_param_x = X, maxIter = maxIter,
                      mG = mG, mE = mE, maxDiag = maxDiag, traits = colnames(Y),
                      tolerance = 1e-4 * max(abs(solve(var(Y)))))
      if (!is.null(snpCov)) {
        ## Without snpCovariates.
        varCompRed <- EMFA(y = Y, k = K, size_param_x = XRed,
                           maxIter = maxIter, mG = mG, mE = mE,
                           maxDiag = maxDiag, traits = colnames(Y))
      }
    } else if (GLSMethod == "multi") {
      ## Including snpCovariates.
      varComp <- sapply(X = chrs, FUN = function(chr) {
        EMFA(y = Y, k = K[[which(chrs == chr)]], size_param_x = X,
             maxIter = maxIter, mG = mG, mE = mE, maxDiag = maxDiag,
             traits = colnames(Y))
      }, simplify = FALSE)
      if (!is.null(snpCov)) {
        ## Without snpCovariates.
        varCompRed <- sapply(X = chrs, FUN = function(chr) {
          EMFA(y = Y, k = K[[which(chrs == chr)]], size_param_x = XRed,
               maxIter = maxIter, mG = mG, mE = mE, maxDiag = maxDiag,
               traits = colnames(Y))
        }, simplify = FALSE)
      }
    }
  }
  if (GLSMethod == "single") {
    res <- list(Vg = varComp$Vg,
                Ve = varComp$Ve,
                VgRed = varCompRed$Vg,
                VeRed = varCompRed$Ve,
                traits = colnames(Y))
  } else if (GLSMethod == "multi") {
    res <- list(Vg = lapply(X = varComp, FUN = `[[`, "Vg"),
                Ve = lapply(X = varComp, FUN = `[[`, "Ve"),
                VgRed = lapply(X = varCompRed, FUN = `[[`, "Vg"),
                VeRed = lapply(X = varCompRed, FUN = `[[`, "Ve"),
                traits = colnames(Y))
  }
  return(res)
}

#' Helper function for restricting traits based on heritability
#'
#' @noRd
#' @keywords internal
restrictTraits <- function(Y,
                           X,
                           K) {
  ## For multiple kinship matrices, combine them.
  if (!is.matrix(K)) {
    K <- Reduce(f = `+`, x = K) / length(K)
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
  nTrait <- ncol(Y)
  ratio <- sapply(X = seq_len(nTrait), FUN = function(i) {
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
                           data = dat, verbose = FALSE, dateWarning = FALSE)
    ## Extract components from fitted model.
    Vg <- as.numeric(modFit$sigma[[1]])
    Ve <- as.numeric(modFit$sigma[[2]])
    return(Vg / (Vg + Ve))
  })
  traitsMod <- traits[ratio > 0.01]
  nTraitMod <- length(traitsMod)
  if (nTraitMod < 2) {
    stop("Only ", nTraitMod, " traits have a heritability larger than 0.01. ",
         "An unstructured model cannot be fitted. Please try a factor ",
         "analytic or pairwise model.\n",
         call. = FALSE)
  }
  if (nTrait != nTraitMod) {
    warning("The following traits have a heritability smaller than 0.01 ",
            "and are removed from the model:\n",
            paste(setdiff(traits, traitsMod), collapse = ", "), "\n",
            call. = FALSE)
  }
  return(traitsMod)
}



#' Select markers to be excluded from GWAS scan.
#'
#' Helper function for selecting markers to be excluded from GWAS scan.
#' Markers are excluded if they are identical to any of the snpCovariates
#' (including the snpCovariates themselves).
#'
#' @param snpCov A character vector of snpCovariates.
#' @param markers A matrix with marker information.
#' @param allFreq A numerical vector of allele frequencies of the markers in
#' \code{markers}. This could be computed from markers as well but it is
#' needed in the general algorithm so to not redo things unnecessarily it is
#' not redone here.
#'
#' @returns A numerical vector of markers to be exluded from the GWAS scan.
#'
#' @keywords internal
exclMarkers <- function(snpCov,
                        markers,
                        allFreq,
                        ref = NULL) {
  exclude <- integer()
  if (any(snpCov %in% colnames(markers))) {
    snpCovNumbers <- which(colnames(markers) %in% snpCov)
    if (length(dim(markers)) == 2) {
      for (snp in snpCovNumbers) {
        ## Rough selection based on allele frequency. Done for speed.
        candidates <- which(allFreq == allFreq[snp])
        ## Exclude all snps that are identical to snps in snpCovariates.
        snpInfo <- as.numeric(markers[, snp])
        exclude <- union(exclude,
                         candidates[apply(X = markers[, candidates,
                                                      drop = FALSE],
                                          MARGIN = 2, FUN = function(x) {
                                            identical(as.numeric(x), snpInfo)
                                          })])
      }
    } else if (length(dim(markers)) == 3) {
      ## Compute mean value for reference allele.
      allMeans <- apply(markers[ , , -ref], c(3, 2), mean)
      for (snp in snpCovNumbers) {
        for (allele in rownames(allMeans[allMeans[, snp] != 0, ])) {
          ## Rough selection based on mean. Done for speed.
          candidates <- which(allMeans == allMeans[allele, snp], arr.ind = TRUE)
          exclude <- union(exclude,
                           candidates[apply(X = candidates, MARGIN = 1,
                                            FUN = function(m) {
                                              identical(markers[, m[2], m[1]],
                                                        markers[, snp, allele])
                                            }), 2])
        }
      }
      ## Rough selection based on mean. Done for speed.
      candidates <- which(allMeans == allMeans[snp])
      ## Exclude all snps that are identical to snps in snpCovariates.
      snpInfo <- as.numeric(markers[, snp, ])
      exclude <- union(exclude,
                       candidates[apply(X = markers[, candidates, ,
                                                    drop = FALSE],
                                        MARGIN = 2, FUN = function(x) {
                                          identical(as.numeric(x), snpInfo)
                                        })])

    }
  }
  return(exclude)
}


