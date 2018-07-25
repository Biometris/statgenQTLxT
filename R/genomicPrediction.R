#' Perform genomic prediction
#'
#' Perform genomic prediction
#'
#' @inheritParams runSingleTraitGwas
#'
#' @param training a character vector with names of genotypes to use as
#' training data for the prediction.
#' @param keep a character vector with names of columns in \code{covar} in
#' \code{gData} that should be included in the output.
#'
#' @return A data.frame containing predictions for genotypes per trait.
#' Columns in \code{keep} are included in the output.
#'
#' @references Piepho, H.-P., Möhring, J., Schulz-Streeck, T. and Ogutu, J. O.
#' (2012), A stage-wise approach for the analysis of multi-environment trials.
#' Biom. J., 54: 844–860. doi:10.1002/bimj.201100219
#'
#' @export
genomicPrediction <- function(gData,
                              traits = NULL,
                              environments = NULL,
                              training,
                              kin = NULL,
                              kinshipMethod = c("astle", "GRM", "IBS",
                                                "vanRaden"),
                              keep = NULL) {
  ## If environments is null set environments to all environments in pheno.
  if (is.null(environments)) {
    environments <- 1:length(gData$pheno)
  }
  if (!is.null(traits) && !is.numeric(traits) && !is.character(traits)) {
    stop("traits should be a numeric or character vector.\n")
  }
  for (environment in environments) {
    if ((is.character(traits) &&
         !all(traits %in% colnames(gData$pheno[[environment]]))) ||
        (is.numeric(traits) &&
         (any(traits == 1) ||
          any(traits > ncol(gData$pheno[[environment]]))))) {
      stop("traits should be columns in pheno.\n")
    }
  }
  if (is.null(gData$kinship) && is.null(kin)) {
    ## Compute kinship matrix.
    kinshipMethod <- match.arg(kinshipMethod)
    KTrain <- do.call(kinshipMethod, list(X = gData$markers[training, ]))
  } else if (is.null(kin)) {
    ## Restrict K to training dataset.
    KTrain <- gData$kinship[training, training]
  } else if (is.matrix(kin)) {
    ## Restrict K to training dataset.
    KTrain <- as(kin, "dsyMatrix")[training, training]
  } else {
    ## Restrict K to training dataset.
    KTrain <- kin[training, training]
  }
  p <- Matrix::colSums(gData$markers) / (2 * nrow(gData$markers))
  Z <- scale(gData$markers, center = 2 * p)
  Z <- Z / (2 * sqrt(sum(p * (1 - p))))
  Z_sel <- Z[training, ]
  ## Extract training data for modeling.
  dataPredict <- gData$pheno[[environments]]
  dataTrain <- dataPredict[dataPredict$genotype %in% training, ]
  ## Calculate predictions per trait.
  preds <- lapply(X = traits, FUN = function(trait) {
    ## Fit mixed model on training data.
    sommerFit <- sommer::mmer2(fixed = as.formula(paste0(trait, "~ 1")),
                               random = ~ g(genotype),
                               G = list(genotype = KTrain),
                               rcov = ~ units,
                               data = dataTrain,
                               silent = TRUE, date.warning = FALSE)
    ## Extract y.
    y <- dataTrain[[trait]]
    ## Extract mu from model.
    mu <- as.numeric(sommerFit$beta.hat)
    ## Extract genetic variance from model.
    sigma2_u <- sommerFit$sigma["g(genotype)"]
    ## Extract inverse of V from model.
    VInv <- sommerFit$V.inv
    ## Compute markereffects.
    u <- sigma2_u * t(Z_sel) %*% (VInv %*% (y - mu))
    ## Compute predictions.
    g <- Z %*% u + mu
    ## Save predictions for the current trait in a data.frame
    pred <- setNames(data.frame(rownames(g), as.vector(g)),
                     c("genotype", trait))
    return(pred)
  })
  ## Create basedata with columns to be copied to output data.
  baseData <- cbind(data.frame(genotype = rownames(gData$covar)),
                    gData$covar[colnames(gData$covar) %in% keep])
  ## Merge predictions to baseData.
  predTot <- Reduce(f = merge, x = preds, init = baseData)
  return(predTot)
}

