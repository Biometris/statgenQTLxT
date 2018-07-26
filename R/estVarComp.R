#' @keywords internal
estVarComp <- function(GLSMethod,
                       remlAlgo,
                       trait,
                       phenoEnvirTrait,
                       covarEnvir,
                       K,
                       chrs,
                       KChr,
                       nonMissing,
                       nonMissingRepId) {
  ## Estimate variance components.
  if (GLSMethod == 1) {
    if (isTRUE(all.equal(K, Matrix::Diagonal(nrow(K)), check.names = FALSE))) {
      ## Kinship matrix is computationally identical to identity matrix.
      vcovMatrix <- Matrix::Diagonal(nrow(phenoEnvirTrait))
    }
  } else if (GLSMethod == 2) {
    varComp <- vcovMatrix <- setNames(vector(mode = "list",
                                             length = length(chrs)),
                                      paste("chr", chrs))
  }
  if (remlAlgo == 1) {
    ## emma algorithm takes covariates from gData.
    gDataEmma <-
      createGData(pheno = phenoEnvirTrait[, c("genotype", trait)],
                  covar = if (is.null(covarEnvir)) {
                    NULL
                  } else {
                    as.data.frame(phenoEnvirTrait[covarEnvir],
                                  row.names = phenoEnvirTrait$genotype)
                  })
    if (GLSMethod == 1) {
      remlObj <- runEmma(gData = gDataEmma, trait = trait, environment = 1,
                         covar = covarEnvir, K = K)
      ## Extract varComp and vcovMatrix
      varComp <- remlObj$varComp
      vcovMatrix <- remlObj$vcovMatrix
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMissing, nonMissing]
        ## Compute variance components using chromosome specific kinship.
        remlObj <- runEmma(gData = gDataEmma, trait = trait,
                           environment = 1, covar = covarEnvir, K = K)
        ## Compute varcov matrix using var components.
        varComp[[which(chrs == chr)]] <- remlObj$varComp
        vcovMatrix[[which(chrs == chr)]] <- remlObj$vcovMatrix
      }
    }
  } else if (remlAlgo == 2) {
    if (!is.null(covarEnvir)) {
      ## Construct the formula for the fixed part of the model.
      ## Define formula for fixed part. ` needed to accommodate -
      ## in variable names.
      fixed <- as.formula(paste0(trait," ~ `", paste0(covarEnvir,
                                                      collapse = "` + `"),
                                 "`"))
    } else {
      fixed <- as.formula(paste(trait, " ~ 1"))
    }
    if (GLSMethod == 1) {
      ## Fit model.
      modFit <- sommer::mmer2(fixed = fixed, data = phenoEnvirTrait,
                              random = ~ g(genotype), G = list(genotype = K),
                              silent = TRUE, date.warning = FALSE)
      ## Compute varcov matrix using var components from model.
      varCompMod <- modFit$var.comp
      sommerK <- K[nonMissingRepId, nonMissingRepId]
      varComp <- setNames(
        unlist(varCompMod)[c(1, length(unlist(varCompMod)))], c("Vg", "Ve"))
      vcovMatrix <- unlist(varCompMod)[1] * sommerK +
        Matrix::Diagonal(n = nrow(sommerK),
                         x = unlist(varCompMod)[length(unlist(varCompMod))])
      if (any(eigen(vcovMatrix, symmetric = TRUE,
                    only.values = TRUE)$values <= 1e-8))
        vcovMatrix <- Matrix::nearPD(vcovMatrix)$mat
    } else if (GLSMethod == 2) {
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        K <- KChr[[which(chrs == chr)]][nonMissing, nonMissing]
        ## Fit mmer2 model using chromosome specific kinship.
        modFit <- sommer::mmer2(fixed = fixed, data = phenoEnvirTrait,
                                random = ~ g(genotype),
                                G = list(genotype = K), silent = TRUE,
                                date.warning = FALSE)
        ## Compute varcov matrix using var components from model.
        varCompMod <- modFit$var.comp
        sommerK <- K[nonMissingRepId, nonMissingRepId]
        varComp[[which(chrs == chr)]] <- setNames(
          unlist(varCompMod)[c(1, length(unlist(varCompMod)))], c("Vg", "Ve"))
        vcovMatrix[[which(chrs == chr)]] <- unlist(varCompMod)[1] * sommerK +
          unlist(varCompMod)[length(unlist(varCompMod))] *
          Matrix::Diagonal(n = nrow(sommerK))
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
