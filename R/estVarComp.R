#' @keywords internal
estVarComp <- function(GLSMethod,
                       remlAlgo,
                       trait,
                       pheno,
                       covar,
                       K,
                       chrs,
                       KChr,
                       nonMiss,
                       nonMissRepId) {
  ## Estimate variance components.
  if (GLSMethod == 1) {
    if (isTRUE(all.equal(K, Matrix::Diagonal(nrow(K)), check.names = FALSE))) {
      ## Kinship matrix is computationally identical to identity matrix.
      vcovMatrix <- Matrix::Diagonal(nrow(pheno))
    }
  } else if (GLSMethod == 2) {
    varComp <- vcovMatrix <-
      setNames(vector(mode = "list", length = length(chrs)), paste("chr", chrs))
  }
  if (remlAlgo == 1) {
    ## emma algorithm takes covariates from gData.
    gDataEmma <-
      createGData(pheno = pheno[, c("genotype", trait)],
                  covar = if (is.null(covar)) {
                    NULL
                  } else {
                    as.data.frame(pheno[covar], row.names = pheno$genotype)
                  })
    if (GLSMethod == 1) {
      remlObj <- runEmma(gData = gDataEmma, trait = trait, environment = 1,
                         covar = covar, K = K)
      ## Extract varComp and vcovMatrix
      varComp <- remlObj$varComp
      vcovMatrix <- remlObj$vcovMatrix
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMiss, nonMiss]
        ## Compute variance components using chromosome specific kinship.
        remlObj <- runEmma(gData = gDataEmma, trait = trait,
                           environment = 1, covar = covar, K = K)
        ## Compute varcov matrix using var components.
        varComp[[which(chrs == chr)]] <- remlObj$varComp
        vcovMatrix[[which(chrs == chr)]] <- remlObj$vcovMatrix
      }
    }
  } else if (remlAlgo == 2) {
    if (!is.null(covar)) {
      ## Construct the formula for the fixed part of the model.
      ## Define formula for fixed part. ` needed to accommodate -
      ## in variable names.
      fixed <- as.formula(paste0(trait," ~ `",
                                 paste0(covar, collapse = "` + `"), "`"))
    } else {
      fixed <- as.formula(paste(trait, " ~ 1"))
    }
    if (GLSMethod == 1) {
      ## Fit model.
      modFit <- sommer::mmer2(fixed = fixed, data = pheno,
                              random = ~ g(genotype), G = list(genotype = K),
                              silent = TRUE, date.warning = FALSE)
      ## Compute varcov matrix using var components from model.
      vcMod <- modFit$var.comp
      modK <- K[nonMissRepId, nonMissRepId]
      varComp <- setNames(
        unlist(vcMod)[c(1, length(unlist(vcMod)))], c("Vg", "Ve"))
      vcovMatrix <- unlist(vcMod)[1] * modK +
        Matrix::Diagonal(n = nrow(modK),
                         x = unlist(vcMod)[length(unlist(vcMod))])
      if (any(eigen(vcovMatrix, symmetric = TRUE,
                    only.values = TRUE)$values <= 1e-8))
        vcovMatrix <- Matrix::nearPD(vcovMatrix)$mat
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMiss, nonMiss]
        ## Fit mmer2 model using chromosome specific kinship.
        modFit <- sommer::mmer2(fixed = fixed, data = pheno,
                                random = ~ g(genotype),
                                G = list(genotype = K), silent = TRUE,
                                date.warning = FALSE)
        ## Compute varcov matrix using var components from model.
        vcMod <- modFit$var.comp
        modK <- K[nonMissRepId, nonMissRepId]
        varComp[[which(chrs == chr)]] <- setNames(
          unlist(vcMod)[c(1, length(unlist(vcMod)))], c("Vg", "Ve"))
        vcovMatrix[[which(chrs == chr)]] <- unlist(vcMod)[1] * modK +
          unlist(vcMod)[length(unlist(vcMod))] *
          Matrix::Diagonal(n = nrow(modK))
      }
      vcovMatrix <- lapply(vcovMatrix, FUN = function(vc) {
        if (any(eigen(vc, symmetric = TRUE,
                      only.values = TRUE)$values <= 1e-8)) {
          Matrix::nearPD(vc)$mat
        } else {
          vc
        }
      })
    }
  }
  return(list(varComp = varComp, vcovMatrix = vcovMatrix))
}
