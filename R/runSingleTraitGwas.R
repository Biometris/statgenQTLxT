#' run single trait GWAS

#' @param gData an object of class \code{gData} containing at least a data.frame \code{pheno}.
#' @param traits a vector of traits on which to run GWAS. These can be either numeric indices or
#' character names of columns in \code{pheno}. If \code{NULL} GWAS is run on all traits.
#' @param fields a vector of fields on which to run GWAS. These can be either numeric indices or
#' character names of list items in \code{pheno}. If \code{NULL} GWAS is run for all fields.
#' @param covar an optional vector of covariates taken into account when running GWAS. These can be either
#' numeric indices or character names of columns in \code{covar} in \code{gData}. If \code{NULL} no
#' covariates are used.
#' @param snpCovariates an optional character vector of snps to be included as covariates.
#' @param remlAlgo an integer indicating the algorithm used to estimate the variance components
#' \enumerate{
#' \item{sommer}
#' \item{emma}
#' }
#' @param GLSMethod an integer indicating the software used to estimate the marker effects.
#' \enumerate{
#' \item{scan_GLS (not available when there are heterozygotes) -- NOT YET IMPLEMENTED}
#' \item{Fast-LMM -- NOT YET IMPLEMENTED}
#' \item{within R}
#' \item{within R, with chromosome specific kinship matrices; similar to the approach of Rincent et al.}
#' }
#' @param sizeIncludedRegion an integer. Should the results for SNPs close to significant SNPs
#' be included? If so, the size of the region in centimorgan or base pairs. Otherwise 0.
#' @param minR2 A numerical value between 0 and 1. Restricts the SNPs included in the region close
#' to significant SNPs to only those SNPs that are in sufficient LD with the significant snp,
#' where LD is measured in terms of r^2. If for example \code{sizeIncludedRegion} = 200000
#' and \code{minR2} = 0.5, then for every significant SNP also those SNPs whose LD (r^2) with the
#' significant SNP is at least 0.5 AND which are at most 200kb away from this significant snp are
#' included. Ignored if \code{sizeIncludedRegion} = 0.
#' @param useMAF should the minor allele frequency be used for selecting SNPs for the analysis.
#' If \code{FALSE} the minor allele count is used.
#' @param MAF A numerical value between 0 and 1. SNPs with minor allele frequency below
#' this value are not taken into account for the analysis, i.e. p-values and effect sizes are
#' put to missing (NA). Ignored if \code{useMAF} is \code{FALSE}.
#' @param MAC A numerical value. SNPs with minor allele count below this value are not
#' taken into account for the analysis, i.e. p-values and effect sizes are put to missing (NA).
#' Ignored if \code{useMAF} is \code{TRUE}.
#' @param genomicControl should genomic control correction as in Devlin and Roeder be applied?
#' @param boundType an integer indicating the type of threshold used for the selection of candidate
#' loci.
#' \enumerate{
#' \item{Bonferroni; a LOD-threshold of -log10(alpha/p), where p is the number of markers and alpha
#' can be specified in alpha.}
#' \item{a self-chosen LOD-threshold, specied in LODThr.}
#' \item{the LOD-threshold is chosen such the the SNPs with the K smallest p-values are selected. K can
#' be specified -- NOT YET IMPLEMENTED}
#' \item{Bonferroni; as option 1 but with p replaced by the number of effective tests approach
#' as in Gao et al. -- NOT YET IMPLEMENTED}
#' }
#' @param alpha a numerical value used for calculating the LOD-threshold for \code{boundType} = 1.
#' @param LODThr a numerical value used as a LOD-threshold when \code{boundType} = 2.
#'
#'
#' @references Rincent et al. (2014) ....
#' @references Devlin, B. and Roeder K. (1999) Genomic control for association studies. Biometrics,
#' December 1999, Vol. 55(4), p. 997-1004.
#' @references Gao et al. (2008, 2010)
#'
#' @import stats

runSingleTraitGwas <- function (gData,
  traits = NULL,
  fields = NULL,
  covar = NULL,
  snpCovariates = NULL,
  remlAlgo = 1,
  GLSMethod = 3,
  sizeIncludedRegion = 0,
  minR2,
  useMAF = TRUE,
  MAF = 0.05,
  MAC = 10,
  genomicControl = FALSE,
  boundType = 1,
  alpha = 0.05 ,
  LODThr = 4) {

  if(!is.null(snpCovariates) && !all(snpCovariates %in% colnames(gData$markers)))
    stop("All snpCovariates should be in markers")

  if (is.null(GLSMethod) || length(GLSMethod) > 1 || !is.numeric(GLSMethod) || !(GLSMethod %in% 1:4)) {
    GLSMethod <- 3
    warning("GLSMethod has to be 1 (scan_GLS), 2 (FasT-LMM), 3 (R) or 4
      (R with chromosome-specific kinship matrices).\n GLSMethod has been set to 3\n")
  }
  if (is.null(sizeIncludedRegion) || length(sizeIncludedRegion) > 1 || !is.numeric(sizeIncludedRegion) ||
      round(sizeIncludedRegion) != sizeIncludedRegion)
    stop("sizeIncludedRegion should be an integer\n")

  if (sizeIncludedRegion > 0) {
    if (is.null(minR2) || length(minR2) > 1 || !is.numeric(minR2))
      stop("minR2 should be an numerical value between 0 and 1\n")
  }

  if (useMAF) {
    if (MAF == 0) {MAF <- 1e-6}
  } else {
    if (MAC == 0) {MAC <- 1}
  }

  if (!remlAlgo %in% c(1, 2)) {remlAlgo <- 2}



  # if (kinship.type==2) {
  #   #genomicControl <- T
  #   gData$external$kinship.name <- alternative.kinship.name
  #   gData$kinship               <- read.table(file=gData$external$kinship.name,sep=",",header=T,check.names=F)
  #   gData$external$kinship.name <- paste("temp_",alternative.kinship.name,sep="")
  #   gData$kinship               <- as.matrix(gData$kinship)
  #   stopifnot(length(setdiff(colnames(gData$kinship),gData$plant.names))==0 & length(setdiff(gData$plant.names,colnames(gData$kinship)))==0)
  #   rownames(gData$kinship)    <- colnames(gData$kinship)
  #   gData$kinship              <- gData$kinship[gData$plant.names,gData$plant.names]
  #   #if (!identical(colnames(gData$kinship),gData$plant.names)) {
  #   #  cat("WARNING: changed column names in imported kinship matrix.","\n")
  #   #  colnames(gData$kinship)   <-  gData$plant.names
  #   #}
  #   gData$kinship.asreml        <- MakeKinshipAsreml(gData$kinship,genotype.names=gData$plant.names)
  # }


  # if (GLSMethod==2) {
  #   if (!file.exists(paste(gData$external$plink.name,".bed",sep=""))) {
  #     cat("Creating PLINK binary files; this may take a while.\n")
  #     MakePlinkFiles(gData,file.name=plink.file.name)
  #     gData <- AddBinaryPlinkFilesToGwasObj(gData=gData,plink.name=plink.file.name)
  #     save(gData,file=r.image.name)
  #     cat("",file="fastlmmLog.txt")
  #     #   if (show.call) {cat(command.string,"\n",file="fastlmmLog.txt")}
  #   }
  #   #!# minor allele fr's should be based on only those accessions for which phenotypic values are available
  #   #minor.allFreqSel <- apply(gData$markers,1,function(x){mean(as.numeric(x))})
  #   #minor.allFreqSel[minor.allFreqSel>.5] <-  1-minor.allFreqSel[minor.allFreqSel>.5]
  #
  # }

  ## Covariates
  field <- 1
  if (is.null(covar)) {
    phenoField <- tibble::rownames_to_column(as.data.frame(gData$pheno[[field]]), var = "genotype")
  } else {
    ## Append covariates to pheno data. Merge to remove values from pheno that are missing in covar.
    phenoField <- merge(gData$pheno[[1]], gData$covar[covar], by = "row.names")
    colnames(phenoField)[1] <- "genotype"
    ## Remove rows from phenoField with missing covar check if there are missing values.
    phenoField <- phenoField[!is.na(phenoField[covar]), ]
    ## Expand covariates that are a factor (i.e. dummy variables are created) using model.matrix
    ## The new dummies are attached to phenoField, and covar is changed accordingly
    factorCovs <- which(sapply(gData$covar[covar], FUN = is.factor))
    if (length(factorCovs) > 0) {
      covFormula <- as.formula(paste("genotype ~ ", paste(covar[factorCovs], collapse = "+")))
      ## Create dummy variables. Remove intercept.
      extraCov <- as.data.frame(suppressWarnings(model.matrix(object = covFormula, data = phenoField))[, -1])
      ## Add dummy variables to pheno data.
      phenoField <- cbind(phenoField[, -which(colnames(phenoField) %in% names(factorCovs))], extraCov)
      ## Modify covar to suit newly defined columns
      covar <- c(covar[-factorCovs], colnames(extraCov))
    }
  }
  if (!is.null(snpCovariates)) {
    ## Add snp covariates to covar.
    covar <- c(covar, snpCovariates)
    ## Add snp covariates to pheno data.
    phenoField <- merge(phenoField, gData$markers[, snpCovariates], by.x = "genotype",
      by.y = "row.names")
    colnames(phenoField)[(ncol(phenoField) - length(snpCovariates) + 1):ncol(phenoField)] <- snpCovariates
  }

  ####################################################

  # if (make.all.trait.significant.snp.file) {cat("",file=paste("results/",all.trait.significant.snp.file,sep=""))}
  #
  # if (make.all.trait.h2.file) {
  #   all.traits <- data.frame(trait=names(gData$pheno)[trait.numbers],h2=rep(NA,length(trait.numbers)),inflation=rep(NA,length(trait.numbers)))
  # }
  #
  # if (make.one.file.with.all.lod.scores) {
  #   all.lod.scores.frame <- data.frame(gData$map,matrix(NA,gData$N,length(trait.numbers)))
  #   names(all.lod.scores.frame)[-(1:ncol(gData$map))] <- names(gData$pheno)[trait.numbers]
  # }
  #
  # if (make.one.file.with.all.snp.effects) {
  #   all.snp.effects.frame <- data.frame(gData$map,matrix(NA,gData$N,length(trait.numbers)))
  #   names(all.snp.effects.frame)[-(1:ncol(gData$map))] <- names(gData$pheno)[trait.numbers]
  # }

  ###################################################################################################


  ## TO DO: loop over fields
  ## TO DO: convert traits to string per field
  ## TO DO: convert covariates to string per field

  ## Compute kinship matrices per chromosome. Only needs to be done once.
  if (GLSMethod == 4) {
    chrs <- unique(gData$map$chr)
    KChr <- vector(mode = "list", length = length(chrs))
    for (chr in chrs) {
      nonChr <- which(colnames(gData$markers) %in% rownames(gData$map[gData$map$chr != chr, ]))
      KChr[[which(chrs == chr)]] <- GRM(gData$markers[, nonChr])
    }
  }
  ## Compute max value in markers
  maxScore <- max(gData$markers, na.rm = TRUE)
  # allele frequencies based on all genotypes (trait-independent)
  allFreqTot <- colMeans(gData$markers, na.rm = TRUE)
  if (maxScore == 2) {
    allFreqTot<- allFreqTot/ 2
  }
  ## If no traits supplied extract them from pheno data. Exclude genotype and covariates
  if (is.null(traits)) traits <- colnames(phenoField[-which(colnames(phenoField) %in% c("genotype", covar))])
  GWATot <- snpSigTot <- NULL
  ## Perform GWAS for all traits.
  for (trait in traits) {
    ## Select relevant columns only.
    phenoFieldTrait <- phenoField[!is.na(phenoField[trait]), c("genotype", trait, covar)]
    ## Select genotypes where trait is not missing.
    nonMissing <- unique(phenoFieldTrait$genotype)
    kinshipRed <- gData$kinship[nonMissing, nonMissing]
    nonMissingRepId <- phenoFieldTrait$genotype
    ## Estimate variance components.
    if (GLSMethod == 3) {
      if (!isTRUE(all.equal(kinshipRed, diag(nrow(kinshipRed)), check.names = FALSE))) {
        if (remlAlgo == 1) {
          ## emma algorithm takes covariates from gData.
          gData$covar <- phenoField[covar]
          remlObj <- runEmma(gData = gData, trait = trait, field = field, covar = covar, K = kinshipRed)
          ## Compute varcov matrix using var components.
          vcovMatrix <- remlObj[[1]][1] * remlObj[[2]] +
            remlObj[[1]][2] * diag(nrow(remlObj[[2]]))
        } else if (remlAlgo == 2) {
          ## Define random part of the model.
          random <- as.formula("~ sommer::g(genotype)")
          ## Construct the formula for the fixed part of the model.
          if (!is.null(covar)) {
            fixed <- as.formula(paste(trait," ~ ", paste(covar, collapse = " + ")))
          } else {
            fixed <- as.formula(paste(trait, " ~ 1"))
          }
          ## Fit mmer2 model.
          sommerFit <- sommer::mmer2(fixed = fixed, data = phenoFieldTrait,
            random = random, G = list(genotype = kinshipRed), silent = TRUE)
          ## Compute varcov matrix using var components from model.
          vcovMatrix <- sommerFit$var.comp[1, 1] * kinshipRed +
            diag(sommerFit$var.comp[nrow(sommerFit$var.comp), 1], nrow = nrow(kinshipRed))
        }
      } else {
        ## Kinship matrix is computationally identical to identity matrix.
        vcovMatrix <- diag(nrow(phenoFieldTrait))
      }
    } else if (GLSMethod == 4) {
      vcovMatrix <- vector(mode = "list", length = length(chrs))
      for (chr in chrs) {
        ## Get chromosome specific kinship.
        KinshipRedChr <- KChr[[which(chrs == chr)]][nonMissing, nonMissing]
        if (remlAlgo == 1) {
          ## emma algorithm takes covariates from gData.
          gData$covar <- phenoField[covar]
          ## Compute variance components using chromosome specific kinship.
          remlObj <- runEmma(gData = gData, trait = trait, field = field,
            covar = covar, K = KinshipRedChr)
          ## Compute varcov matrix using var components.
          vcovMatrix[[which(chrs == chr)]] <- remlObj[[1]][1] * remlObj[[2]] +
            remlObj[[1]][2] * diag(nrow(remlObj[[2]]))
        } else if (remlAlgo == 2) {
          ## Define random part of the model.
          random <- as.formula("~ sommer::g(genotype)")
          ## Construct the formula for the fixed part of the model.
          if (!is.null(covar)) {
            fixed <- as.formula(paste(trait," ~ "), paste(covar, collapse = " + "))
          } else {
            fixed <- as.formula(paste(trait, " ~ 1"))
          }
          ## Fit mmer2 model using chromosome specific kinship.
          sommerFit <- sommer::mmer2(fixed = fixed, data = phenoFieldTrait,
            random = random, G = list(genotype = KinshipRedChr), silent = TRUE)
          ## Compute varcov matrix using var components from model.
          vcovMatrix[[which(chrs == chr)]] <- sommerFit$var.comp[1, 1] * KinshipRedChr +
            diag(sommerFit$var.comp[nrow(sommerFit$var.comp), 1], nrow = nrow(KinshipRedChr))
        }
      }
    }
    ## Compute allele frequencies based on genotypes for which phenotypic data is available.
    markersRed <- gData$markers[nonMissing, ]
    allFreq <- colMeans(markersRed)
    if (maxScore == 2) {
      allFreq <- allFreq / 2
    }
    if (!useMAF) {MAF <- MAC / length(nonMissing) - 1e-5}
    ## Determine segregating markers. Exclude snps used as covariates.
    segMarkers <- which(allFreq >= maxScore * MAF & allFreq <= maxScore * (1 - MAF))
    ## Exclude snpCovariates from segragating markers.
    if (!is.null(snpCovariates)) {
      snpCovariateNumbers <- which(colnames(markersRed) %in% snpCovariates)
      for (snp in snpCovariateNumbers) {
        ## Rough selection based on allele frequency. Done for speed.
        candidates <- which(allFreq == allFreq[snp])
        ## Exclude all snps that are identical to snps in snpCovariates.
        snpInfo <- markersRed[, snp]
        exclude <- candidates[apply(markersRed[, candidates], 2,
          function(x) {identical(as.numeric(x), as.numeric(snpInfo))})]
      }
      ## Remove excluded snps from segregating markers.
      segMarkers <- setdiff(segMarkers, exclude)
    }
    ## Define data.frame for results.
    GWAResult <- data.frame(snp = rownames(gData$map),
      gData$map,
      pValue = NA,
      effect = NA,
      effectSe = NA,
      RLR2 = NA,
      allFreq = allFreq,
      row.names = rownames(gData$map),
      stringsAsFactors = FALSE)
    y <- phenoFieldTrait[, trait]
    if (GLSMethod == 3) {
      ## The following is based on the genotypes, not the replicates:
      X <- markersRed[nonMissingRepId, segMarkers]
      if (length(covar) == 0) {
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix)
      } else {
        ## Define covariate matrix Z.
        Z <- as.matrix(phenoFieldTrait[which(phenoFieldTrait$genotype %in% nonMissing), covar])
        ## Compute pvalues and effects using fastGLS.
        GLSResult <- fastGLSCov(y = y, X = X, Sigma = vcovMatrix, covs = Z)
      }
      GWAResult[segMarkers, c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult
    } else if (GLSMethod == 4) {
      ## Same as for GLSMethod 3 but then per chromosome using the chromosome specific
      ## kinship calculated before.
      for (chr in chrs) {
        segMarkersChr <- intersect(segMarkers, which(gData$map$chr == chr))
        X <- markersRed[nonMissingRepId, segMarkersChr]
        if (length(covar) == 0) {
          GLSResult <- fastGLS(y = y, X = X, Sigma = vcovMatrix[[which(chrs == chr)]])
        } else {
          Z <- as.matrix(phenoFieldTrait[which(phenoFieldTrait$genotype %in% nonMissing), covar])
          GLSResult <- fastGLSCov(y = y, X = X, Sigma = vcovMatrix[[which(chrs == chr)]], covs = Z)
        }
        GWAResult[segMarkersChr, c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult
      }
    }
    ## Effects should be for a single allele, not for 2
    if (maxScore == 1) {
      GWAResult$effect <- 0.5 * GWAResult$effect
    }
    ## Calculate the genomic inflation factor and rescale p-values
    if (genomicControl) {
      GC <- genomicControlPValues(pVals = GWAResult$pValue,
        nObs = length(nonMissing),
        nCov = length(covar))
      GWAResult$pValue <- GC[[1]]
    }
    ## Add gene information if available.
    if ("genes" %in% names(gData)) {
      GWAResult <- cbind(GWAResult, gene1 = gData$map$gene1, gene2 = gData$map$gene2)
    }
    ## When boundType is 1, 3 or 4, determine the LOD threshold.
    if (boundType == 1) {
      ## Compute LOD using Bonferroni correction.
      nEff <- sum(!is.na(GWAResult$pValue))
      LODThr <- -log10(alpha / nEff)
      # } else if (boundType == 3) {
      #     LODThr <- sort(-log10(na.omit(GWAResult$pvalue)),decreasing=T)[K]
      # } else if (boundType == 4) {
      #   cut.off <- 0.995
      #   number.of.nonmissing  <- aggregate(gData$pheno[,trait],by=list(ordered(gData$pheno$genotype)),FUN=function(x){sum(!is.na(x))})[match(gData$plant.names,sort(gData$plant.names)),2]
      #   number.of.nonmissing[number.of.nonmissing>0] <- 1
      #   ind.indices <- rep((1:gData$n)[number.of.nonmissing>0],times=number.of.nonmissing[number.of.nonmissing>0])
      #   SIGMA       <- varcomp.values[1,] * gData$kinship[ind.indices,ind.indices] + varcomp.values[2,] * diag(sum(number.of.nonmissing))
      #   COR         <- cov2cor(SIGMA)
      #   INV.COR     <- GINV(COR)
      #   Keff      <- 0
      #   n.block   <- 0
      #   b.size    <- 10*sum(number.of.nonmissing)
      #   for (CHR in 1:gData$nchr) {
      #     blocks  <- DefineBlocks(which(gData$map$chromosome==CHR),block.size=b.size)
      #     n.block <- n.block  +  length(blocks)
      #     for (b in 1:length(blocks)) {
      #       marker.frame <- gData$markers[blocks[[b]],]
      #       Keff <- Keff + GaoCorrection(marker.frame=marker.frame,number.of.replicates=number.of.nonmissing,inv.cor.matrix=INV.COR,cut.off=cut.off)
      #     }
      #   }
      #   LODThr <- -log10(alpha/Keff)
    }

    ##Select the SNPs whose LOD-scores is above the threshold
    signSnp <- which(!is.na(GWAResult$pValue) & -log10(GWAResult$pValue) >= LODThr)
    if (length(signSnp) > 0) {
      if (sizeIncludedRegion > 0) {
        snpSelection <- unlist(sapply(signSnp,
          FUN = getSNPsInRegionSufLD,
          gData = gData,
          regionSize = sizeIncludedRegion,
          minR2 = minR2))
        snpSelection <- sort(union(snpSelection, signSnp))
        snpStatus <- rep(paste("within", sizeIncludedRegion/1000 , "kb of a significant snp"),
          length(snpSelection))
        snpStatus[snpSelection %in% signSnp] <- "significant snp"
      } else {
        snpSelection <- signSnp
        snpStatus <- rep("significant snp", length(signSnp))
      }
      if (GLSMethod %in% 3:4) {
        effects <- GWAResult$effect[snpSelection]
        ## Compute variance of marker scores, based on genotypes for which phenotypic data is available.
        ## for inbreeders, this depends on maxScore. It is therefore scaled to marker scores 0, 1
        ## (or 0, 0.5, 1 if there are heterozygotes)
        snpVar <- 4 * effects ^ 2 * apply(markersRed[, snpSelection], 2, var) / maxScore ^ 2
        propSnpVar <- snpVar / as.numeric(var(phenoFieldTrait[trait]))
      }
      ## Create data.frame with significant snps.
      snpSig <- data.frame(marker = colnames(gData$markers)[snpSelection],
        gData$map[snpSelection, ],
        pValue = GWAResult$pValue[snpSelection],
        snpStatus,
        allFreq = allFreq[snpSelection],
        effects,
        propVarLRT = GWAResult$RLR2[snpSelection],
        propSnpVar = propSnpVar,
        stringsAsFactors = FALSE)
      snpSigTot <- rbind(snpSigTot, data.frame(trait = trait, snpSig, stringsAsFactors = FALSE))
    }
    GWATot <- rbind(GWATot, data.frame(trait = trait, GWAResult, stringsAsFactors = FALSE))


    ############################################################################

    # Create the summary file:
    #
    #     summ.file    <- paste("results/","summary.",trait,suffix,".txt",sep="")
    #
    #
    #     cat("R-image used: ",r.image.name,"\n","\n",file=summ.file)     # file is created and first line is written; APPEND=FALSE (default)
    #     cat("Trait: ",trait,"\n","\n",file=summ.file,append=T)
    #     cat("Analysis started on: ",start.date,"\n",file=summ.file,append=TRUE)
    #     cat("Analysis finished on: ",date(),"\n","\n",file=summ.file,append=TRUE)
    #
    #     cat("Data are available for",gData$N,"SNPs", "\n",file=summ.file,append=TRUE)
    #     if (MAF > 0) {cat(gData$N-nEff,"of them were not analyzed because their minor allele frequency is below",MAF,"\n","\n",file=summ.file,append=TRUE)}
    #
    #
    #     if (GLSMethod %in% c(1,3)) {
    #       cat("Mixed model with only polygenic effects, and no marker effects:","\n",file=summ.file,append=TRUE)
    #       cat("Genetic variance: ",varcomp.values[1,1],"\t","standard error: ",varcomp.std[1],"\n",file=summ.file,append=TRUE)
    #       cat("Residual variance: ",varcomp.values[nrow(varcomp.values),1],"\t","standard error: ",varcomp.std[nrow(varcomp.values)],"\n\n",file=summ.file,append=TRUE)
    #       #
    #       cat("File containing the p-values of all snps:  ",output.file,"\n","\n",append=TRUE,file=summ.file,sep="")
    #     }
    #     #
    #
    #     if (boundType %in% 1:4) {
    #       cat("LOD-threshold: ",LODThr,"\n",append=TRUE,file=summ.file)
    #       if (boundType==4) {cat("Number of effective tests: ",Keff,"\n",append=TRUE,file=summ.file)}
    #       #
    #       if (!no.significant.snps) {
    #         cat("File containing the p-values of the selected snps: ",snp.output.file,"\n",append=TRUE,file=summ.file)
    #         cat("Number of selected snps =",nrow(snp.sig),"\n",append=TRUE,file=summ.file)
    #         cat("Smallest p-value among the selected snps:",min(snp.sig$pvalue),"\n",append=TRUE,file=summ.file)
    #         cat("Largest  p-value among the selected snps:",max(snp.sig$pvalue),"(LOD-score:",-log10(max(snp.sig$pvalue)),")","\n",append=TRUE,file=summ.file)
    #       } else {
    #         cat("No significant snps found.","\n",append=TRUE,file=summ.file)
    #       }
    #       if (genomicControl) {cat("\n","Genomic control correction was applied","\n",append=TRUE,file=summ.file)}
    #       if (!genomicControl) {cat("\n","No Genomic control correction","\n",append=TRUE,file=summ.file)}
    #       cat("Genomic control inflation-factor = ",GC[[2]],"\n","\n",append=TRUE,file=summ.file)
    #     }




  } # end for (trait in traits)
  return(list(GWAResult = GWATot, signSnp = snpSigTot))
}
