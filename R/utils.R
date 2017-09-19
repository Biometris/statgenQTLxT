#' @keywords internal

### Add covariates and snpCovariates to phenotypic data and convert covariate factors to dummy varables.
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
    phenoEnvir <- merge(phenoEnvir, as.matrix(gData$markers[, snpCovariates, drop = FALSE]),
      by.x = "genotype", by.y = "row.names")
    colnames(phenoEnvir)[(ncol(phenoEnvir) - length(snpCovariates) + 1):ncol(phenoEnvir)] <- snpCovariates
  }
  return(list(phenoEnvir = phenoEnvir, covarEnvir = covarEnvir))
}

## Compute chromosome specific kinship matrices.
chrSpecKin <- function(gData, kinshipMethod) {
  chrs <- unique(gData$map$chr[rownames(gData$map) %in% colnames(gData$markers)])
  if (length(chrs) == 1)
    stop("Chromosome specific kinship calculation not possible since map contains
      only 1 chromosome.\n")
  ## Create list of zero matrices.
  KChr <- setNames(replicate(n = length(chrs),
    Matrix::Matrix(data = 0, nrow = nrow(gData$markers), ncol = nrow(gData$markers),
      dimnames = list(rownames(gData$markers), rownames(gData$markers))),
    simplify = FALSE),
    paste0("KChr", chrs))
  ## Create vector of marker numbers per chromosome.
  nMrkChr <- setNames(numeric(length = length(chrs)), chrs)
  for (chr in chrs) {
    ## Extract markers for current chromosome.
    chrMrk <- which(colnames(gData$markers) %in% rownames(gData$map[gData$map$chr == chr, ]))
    ## Compute kinship for current chromosome only. Denominator = 1, division is done later.
    K <- do.call(kinshipMethod, list(X = gData$markers[, chrMrk, drop = FALSE], denominator = 1))
    ## Compute number of markers for other chromosomes.
    nMrkChr[which(chrs == chr)] <- ncol(gData$markers[, -chrMrk, drop = FALSE])
    ## Add computed kinship to all other matrices in KChr.
    for (i in setdiff(1:length(chrs), which(chr == chrs))) {
      KChr[[i]] <- KChr[[i]] + K
    }
  }
  ## Divide matrix for current chromosome by number of markers in other chromosomes.
  for (i in 1:length(KChr)) {
    KChr[[i]] <- KChr[[i]] / nMrkChr[i]
  }
  return(KChr)
}

## Select markers to be excluded from GWAS scan.
computeExcludedMarkers <- function(snpCovariates, markersRed, allFreq) {
  exclude <- integer()
  if (any(snpCovariates %in% colnames(markersRed))) {
    snpCovariateNumbers <- which(colnames(markersRed) %in% snpCovariates)
    for (snp in snpCovariateNumbers) {
      ## Rough selection based on allele frequency. Done for speed.
      candidates <- which(allFreq == allFreq[snp])
      ## Exclude all snps that are identical to snps in snpCovariates.
      snpInfo <- markersRed[, snp]
      exclude <- union(exclude, candidates[apply(markersRed[, candidates], MARGIN = 2,
        FUN = function(x) {identical(as.numeric(x), as.numeric(snpInfo))})])
    }
  }
  return(exclude)
}

## Fill GWAResult data.frame for (a selection of) markers
fillGWAResult <- function(GWAResult, effects, effectsSe, Xt, Yt, VInvArray,
  excludedMarkers, markersRed, Uk) {
  p <- ncol(effects)
  est0 <- estimateEffects(X = Xt, Y = Yt, VInvArray = VInvArray, returnAllEffects = TRUE)
  fittedMean0 <- Matrix::Matrix(est0$effectsEstimates, ncol = length(est0$effectsEstimates) / p) %*% Xt
  SS0 <- LLQuadFormDiag(Y = Yt - fittedMean0, VInvArray = VInvArray)
  for (mrk in setdiff(1:ncol(markersRed), excludedMarkers)) {
    mrkName <- colnames(markersRed)[mrk]
    x <- markersRed[, mrk, drop = FALSE]
    xt <- Matrix::crossprod(x, Uk)
    LRTRes <- LRTTest(X = Xt, x = xt, Y = Yt, VInvArray = VInvArray, SS0 = SS0)
    GWAResult[mrkName, "pValue"] <- LRTRes$pvalue
    GWAResult[mrkName, "pValueWald"] <- pchisq(sum((LRTRes$effects / LRTRes$effectsSe) ^ 2),
      df = p, lower.tail = FALSE)
    effects[mrkName, ] <- LRTRes$effects
    effectsSe[mrkName, ] <-  LRTRes$effectsSe
  }
  return(list(GWAResult = GWAResult, effects = effects, effectsSe = effectsSe))
}





