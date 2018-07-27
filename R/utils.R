## Add covariates and snpCovariates to phenotypic data and convert covariate
## factors to dummy varables.
#' @keywords internal
expandPheno <- function(gData,
                        env,
                        covar,
                        snpCov) {
  ## Add covariates to pheno data.
  if (is.null(covar)) {
    phEnv <- gData$pheno[[env]]
    covEnv <- NULL
  } else {
    ## Append covariates to pheno data. Merge to remove values from pheno that
    ## are missing in covar.
    phEnv <- merge(gData$pheno[[env]], gData$covar[covar],
                   by.x = "genotype", by.y = "row.names")
    ## Remove rows from phEnv with missing covar check if there are
    ## missing values.
    phEnv <- phEnv[complete.cases(phEnv[covar]), ]
    ## Expand covariates that are a factor (i.e. dummy variables are created)
    ## using model.matrix. The new dummies are attached to phEnv, and covar
    ## is changed accordingly
    factorCovs <- which(sapply(X = gData$covar[covar], FUN = is.factor))
    if (length(factorCovs) > 0) {
      ## Create dummy variables without intercept.
      covFormula <- as.formula(paste("genotype ~ ",
                                     paste(covar[factorCovs], collapse = "+")))
      extraCov <- as.data.frame(suppressWarnings(
        model.matrix(object = covFormula, data = droplevels(phEnv))))[, -1]
      ## Add dummy variables to pheno data.
      phEnv <- cbind(phEnv[, -which(colnames(phEnv) %in% names(factorCovs))],
                     extraCov)
      ## Modify covar to suit newly defined columns
      covEnv <- c(covar[-factorCovs], colnames(extraCov))
    } else {
      covEnv <- covar
    }
  }
  if (!is.null(snpCov)) {
    ## Add snp covariates to covar.
    covEnv <- c(covEnv, snpCov)
    ## Add snp covariates to pheno data.
    phEnv <- merge(phEnv, as.matrix(gData$markers[, snpCov, drop = FALSE]),
                   by.x = "genotype", by.y = "row.names")
    colnames(phEnv)[(ncol(phEnv) - length(snpCov) + 1):ncol(phEnv)] <- snpCov
  }
  return(list(phEnv = phEnv, covEnv = covEnv))
}

## Helper function for computing (or extracting kinship matrices)
## 1 - If kin is supplied use kin
## 2 - Get kin from gData object
## 3 - Compute kin from markers (and map for GLSMethod 2)
computeKin <- function(GLSMethod,
                       kin,
                       gData,
                       markers,
                       map,
                       kinshipMethod) {
  if (GLSMethod == 1) {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to dsyMatrix.
      K <- as(kin, "dsyMatrix")
    } else {
      if (!is.null(gData$kinship)) {
        ## Get kin from gData object.
        K <- gData$kinship
      } else {
        ## Compute K from markers.
        K <- do.call(kinshipMethod, list(X = markers))
      }
    }
  } else if (GLSMethod == 2) {
    if (!is.null(kin)) {
      ## kin is supplied as input. Convert to dsyMatrices.
      K <- lapply(X = kin, FUN = as, Class = "dsyMatrix")
    } else {
      ## Get kin from gData object.
      if (!is.null(gData$kinship)) {
        K <- gData$kinship
      } else {
        ## Compute chromosome specific kinship matrices.
        K <- chrSpecKin(gData = createGData(geno = markers, map = map),
                        kinshipMethod = kinshipMethod)
      }
    }
  }
  return(K)
}

## Compute chromosome specific kinship matrices.
chrSpecKin <- function(gData,
                       kinshipMethod) {
  chrs <- unique(gData$map$chr[rownames(gData$map) %in%
                                 colnames(gData$markers)])
  if (length(chrs) == 1) {
    stop(paste("Chromosome specific kinship calculation not possible since",
               "map contains only 1 chromosome.\n"))
  }
  ## Create list of zero matrices.
  KChr <- setNames(
    replicate(n = length(chrs),
              Matrix::Matrix(data = 0, nrow = nrow(gData$markers),
                             ncol = nrow(gData$markers),
                             dimnames = list(rownames(gData$markers),
                                             rownames(gData$markers))),
              simplify = FALSE),
    paste0("KChr", chrs))
  ## Create vector of marker numbers per chromosome.
  nMrkChr <- setNames(numeric(length = length(chrs)), chrs)
  for (chr in chrs) {
    ## Extract markers for current chromosome.
    chrMrk <- which(colnames(gData$markers) %in%
                      rownames(gData$map[gData$map$chr == chr, ]))
    ## Compute kinship for current chromosome only. Denominator = 1, division
    ## is done later.
    K <- do.call(kinshipMethod, list(X = gData$markers[, chrMrk, drop = FALSE],
                                     denominator = 1))
    ## Compute number of markers for other chromosomes.
    nMrkChr[which(chrs == chr)] <- ncol(gData$markers[, -chrMrk, drop = FALSE])
    ## Add computed kinship to all other matrices in KChr.
    for (i in setdiff(1:length(chrs), which(chr == chrs))) {
      KChr[[i]] <- KChr[[i]] + K
    }
  }
  ## Divide matrix for current chromosome by number of markers in other
  ## chromosomes.
  for (i in 1:length(KChr)) {
    KChr[[i]] <- KChr[[i]] / nMrkChr[i]
  }
  return(KChr)
}

## Select markers to be excluded from GWAS scan.
computeExcludedMarkers <- function(snpCov,
                                   markersRed,
                                   allFreq) {
  exclude <- integer()
  if (any(snpCov %in% colnames(markersRed))) {
    snpCovNumbers <- which(colnames(markersRed) %in% snpCov)
    for (snp in snpCovNumbers) {
      ## Rough selection based on allele frequency. Done for speed.
      candidates <- which(allFreq == allFreq[snp])
      ## Exclude all snps that are identical to snps in snpCovariates.
      snpInfo <- as.numeric(markersRed[, snp])
      exclude <- union(exclude,
                       candidates[apply(X = markersRed[, candidates],
                                        MARGIN = 2, FUN = function(x) {
                                          identical(as.numeric(x), snpInfo)
                                        })])
    }
  }
  return(exclude)
}

## Helper function for accessing parallel computing functions.
getOper <- function(x) {
  if (x) {
    `%dopar%`
  } else {
    `%do%`
  }
}
