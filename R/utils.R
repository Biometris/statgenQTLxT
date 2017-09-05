#' @keywords internal
expandPheno <- function(gData,
  environment,
  covar,
  snpCovariates) {
  ## Add covariates to pheno data.
  if (is.null(covar)) {
    phenoEnvir <- gData$pheno[[environment]]
    covarEnvir <- NULL
  } else {
    ## Append covariates to pheno data. Merge to remove values from pheno that are missing in covar.
    phenoEnvir <- merge(gData$pheno[[environment]], gData$covar[covar],
      by.x = "genotype", by.y = "row.names")
    ## Remove rows from phenoEnvir with missing covar check if there are missing values.
    phenoEnvir <- phenoEnvir[complete.cases(phenoEnvir[covar]), ]
    ## Expand covariates that are a factor (i.e. dummy variables are created) using model.matrix
    ## The new dummies are attached to phenoEnvir, and covar is changed accordingly
    factorCovs <- which(sapply(X = gData$covar[covar], FUN = is.factor))
    if (length(factorCovs) > 0) {
      ## Create dummy variables without intercept.
      covFormula <- as.formula(paste("genotype ~ ", paste(covar[factorCovs], collapse = "+")))
      extraCov <- as.data.frame(suppressWarnings(model.matrix(object = covFormula,
        data = droplevels(phenoEnvir))))[, -1]
      ## Add dummy variables to pheno data.
      phenoEnvir <- cbind(phenoEnvir[, -which(colnames(phenoEnvir) %in% names(factorCovs))], extraCov)
      ## Modify covar to suit newly defined columns
      covarEnvir <- c(covar[-factorCovs], colnames(extraCov))
    } else {
      covarEnvir <- covar
    }
  }
  if (!is.null(snpCovariates)) {
    ## Add snp covariates to covar.
    covarEnvir <- c(covarEnvir, snpCovariates)
    ## Add snp covariates to pheno data.
    phenoEnvir <- merge(phenoEnvir, gData$markers[, snpCovariates], by.x = "genotype",
      by.y = "row.names")
    colnames(phenoEnvir)[(ncol(phenoEnvir) - length(snpCovariates) + 1):ncol(phenoEnvir)] <- snpCovariates
  }
  return(list(phenoEnvir = phenoEnvir, covarEnvir = covarEnvir))
}
