#' run single trait GWAS

#' @param gData an object of class \code{gData} containing at least a data.frame \code{pheno}.
#' @param traits a vector of traits on which to run GWAS. These can be either numeric indices or
#' character names of columns in \code{pheno}. If \code{NULL} GWAS is run on all traits.
#' @param fields a vector of fields on which to run GWAS. These can be either numeric indices or
#' character names of list items in \code{pheno}. If \code{NULL} GWAS is run for all fields.
#' @param covar an optional vector of covariates taken into account when running GWAS. These can be either
#' numeric indices or character names of columns in \code{covar} in \code{gData}. If \code{NULL} no
#' covariates are used.
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

  # Create the subdirectory structure. If it already exists, nothing is done.
  #dir.create("results", showWarnings = F)
  #dir.create("plots", showWarnings = F)
  #dir.create("output", showWarnings = F)

  ####################################################################################################################################################
  # load functions and data:

  if (GLSMethod == 4) {
    if (!(all(paste0('kchr',gData$chromosomes) %in% names(gData)))) {
      stop("If GLSMethod=4, gData should contain chromosome-specific-kinship matrices kchr1, kchr2, etc.")
    }
    if (kinship.type==2) {
      kinship.type <- 1
      cat(paste0('Warning: GLSMethod=4: Therefore the kinship matrix contained in ',alternative.kinship.name,' is ignored.','\n'))
    }

    # if (make.all.trait.h2.file==TRUE) {
    #    make.all.trait.h2.file <- FALSE
    #    cat(paste0('Warning: GLSMethod=4: Therefore make.all.trait.h2.file is set to FALSE','\n'))
    #  }
  }


  if (!all(trait.numbers %in% setdiff(2:ncol(gData$pheno),which.columns.as.factor))) {
    stop("ERROR: all numbers in the vector trait.numbers should be between 2 and the \n number of columns in the .csv file with the phenotypes. \n In addition, colmun-numbers declared in which.columns.as.factor cannot be in trait.numbers")
  }

  ##################################################################################################################

  if (!remlAlgo %in% c(1, 2)) {remlAlgo <- 2}

  # For all columns in gData$pheno that are factors, check if there are missing values; in that case, remove the corresponding rows
  if (covariables) {

    cov.cols.is.factor <- (lapply(as.data.frame(gData$pheno[,cov.cols]),FUN=class)=="factor")

    if (sum(cov.cols.is.factor)>0) {
      for (j in cov.cols[cov.cols.is.factor]) {

        gData$pheno <- gData$pheno[!is.na(gData$pheno[,j]),]

      }
    }

    # also remove rows with missing values for covariates that are no a factor
    if ((length(cov.cols) - sum(cov.cols.is.factor)) > 0) {
      for (j in cov.cols[!cov.cols.is.factor]) {
        gData$pheno <- gData$pheno[!is.na(gData$pheno[,j]),]
      }
    }

  }

  if (kinship.type==2) {
    #genomicControl <- T
    gData$external$kinship.name <- alternative.kinship.name
    gData$kinship               <- read.table(file=gData$external$kinship.name,sep=",",header=T,check.names=F)
    gData$external$kinship.name <- paste("temp_",alternative.kinship.name,sep="")
    gData$kinship               <- as.matrix(gData$kinship)
    stopifnot(length(setdiff(colnames(gData$kinship),gData$plant.names))==0 & length(setdiff(gData$plant.names,colnames(gData$kinship)))==0)
    rownames(gData$kinship)    <- colnames(gData$kinship)
    gData$kinship              <- gData$kinship[gData$plant.names,gData$plant.names]
    #if (!identical(colnames(gData$kinship),gData$plant.names)) {
    #  cat("WARNING: changed column names in imported kinship matrix.","\n")
    #  colnames(gData$kinship)   <-  gData$plant.names
    #}
    gData$kinship.asreml        <- MakeKinshipAsreml(gData$kinship,genotype.names=gData$plant.names)
  }


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

  ################################################################
  ####### covariates

  # if no covariates are to be used, cov.frame will be an dataframe with one column ("mu": vector of ones).
  # Otherwise, it will be a dataframe with ones in the first
  # column, and the actual covariates in subsequent columns

  if (!is.null(covar)) {
    cov.cols  <- integer(0)
    covariate.file  <- ""
  } else {
    # IF ONE OR MORE OF THE COVARIATES IN cov.cols is a factor: these are now "expanded" (i.e. dummy variables are created) using model.matrix
    # The new dummies are attached to gData$pheno, and cov.cols is changed accordingly
    factor.covs <- as.integer(which(lapply(gData$pheno,class)[cov.cols]=="factor"))
    if (sum(factor.covs) >0) {
      temp.formula        <- as.formula(paste(paste("genotype~",paste(names(gData$pheno)[cov.cols[factor.covs]],collapse="+"))))
      extra.cov.frame     <- as.data.frame(suppressWarnings(model.matrix(object=temp.formula,data=gData$pheno))[,-1])
      n.col.old           <- ncol(gData$pheno)
      gData$pheno      <- cbind(gData$pheno,extra.cov.frame)
      #names(gData$pheno)[-(1:n.col.old)] <- names(gData$pheno)[cov.cols[factor.covs]]
      cov.cols            <- c(cov.cols,(ncol(gData$pheno)-ncol(extra.cov.frame)+1):ncol(gData$pheno))
      cov.cols            <- setdiff(cov.cols,cov.cols[factor.covs])
    }
    covariate.file  <- "output/covariates.txt"
  }
  cov.frame       <- MakeCovariateFile(pheno.dataframe=gData$pheno,cov.cols=cov.cols,file.name=covariate.file)


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


  if (!is.null(covar)) {
    temp.pheno <- gData$pheno
  }

  ## TO DO: loop over fields
  ## TO DO: convert traits to string per field
  ## TO DO: convert covariates to string per field

  field <- 1
  phenoField <- tibble::rownames_to_column(as.data.frame(gData$pheno[[field]]), var = "genotype")

  ## Perform GWAS for all traits
  for (trait in traits) {
    # trait <- "anthesis.ARIHAS_2013_drought"

    start.date <- date()

    if (!is.null(covar)) {
      ## if covariates are used, all phenotypic values (for the current trait) for which at least one
      ## covariate is missing, are set to missing
      phenoField[apply(as.data.frame(phenoField[cov.cols]), 1, FUN = anyNA), trait]  <- NA
    }

    ## Select genotypes where trait is not missing
    nonMissing <- unique(phenoField$genotype[!is.na(phenoField[trait])])
    kinshipRed <- gData$kinship[nonMissing, nonMissing]
    nonMissingRepId <- phenoField$genotype[!is.na(phenoField[trait])]

    ## Estimate variance components
    if (remlAlgo == 2) {
      ## Define random part of the model.
      random <- as.formula("~ sommer::g(genotype)")
      ## Construct the formula for the fixed part of the model.
      if (!is.null(covar)) {
        fixed <- as.formula(paste(trait," ~ "), paste(covar, collapse = " + "))
      } else {
        fixed <- as.formula(paste(trait, " ~ 1"))
      }
      ## Fit mmer2 model
      sommerFit <- sommer::mmer2(fixed = fixed, data = phenoField[!is.na(phenoField[, trait]), ],
        random = random, G = list(genotype = kinshipRed), silent = TRUE)
      ## Extract variance components from fitted model.
      varComp <- sommerFit$var.comp[1]
    } else if (remlAlgo == 1 & GLSMethod != 4) {
      emma.obj <- runEmma(gData = gData, trait = trait, field = field, covar = covar, K = K)
      varComp <- data.frame(var.comp.values = emma.obj[[1]])
    }

    #!# to do : if the kinship-matrix is the identity matrix ...


    if (remlAlgo=="emma" & GLSMethod==4) {
      emma.obj.list       <- list()
      varcomp.values.list <- list()

      for (chr in gData$chromosomes) {

        emma.obj.list[[which(chr==gData$chromosomes)]]  <- RunEmma(gwas.obj = gData,trait,cov.cols,
          K.user=gData[[which(names(gData)==paste0('kchr',chr))]])
        varcomp.values.list[[which(chr==gData$chromosomes)]]      <- data.frame(var.comp.values=emma.obj.list[[which(chr==gData$chromosomes)]][[1]])
      }
    }


    ##########################################################
    # CREATE PHENO-TYPE FILE + VARCOMP-FILE FOR  SCAN-GLS / FAST-LMM
    #

    # if (GLSMethod==1) {
    #
    #   varcomp.file        <- paste("output/",trait,".","varcomp",".csv",sep="")
    #
    #   # Name of the phenotype-file for scan_GLS:
    #   input.pheno <- paste("output/",trait,".csv",sep="")
    #
    #   if (sum(abs(gData$kinship-diag(nrow(gData$kinship)))) > 0.0001 ) { # if the kinship-matrix is (numerically) different from the identity matrix ...
    #     # to do : the case where the kinship matrix is diagonal, but has unequal variances
    #     vcovMatrix <- MakeScanGlsKinship(varcomp.values[1,],varcomp.values[nrow(varcomp.values),],
    #       gData$kinship,plant.names=gData$plant.names,
    #       gData$pheno,tr.n = trait)
    #   } else {  # if the kinship-matrix is the identity matrix ...
    #     vcovMatrix <- diag(sum(!is.na(gData$pheno[, trait])))
    #   }
    #   write.table(vcovMatrix,file=varcomp.file,sep=",",quote=F,col.names=F,row.names=F)
    #   MakePhenoFile(pheno.object=gData$pheno,col.number=trait,file.name=input.pheno)
    # }

    if (GLSMethod == 3) {
      if (sum(abs(gData$kinship - diag(nrow(gData$kinship)))) > 0.0001) {
        ## if the kinship-matrix is (numerically) different from the identity matrix ...
        ## TO DO: the case where the kinship matrix is diagonal, but has unequal variances
        vcovMatrix <- emma.obj[[1]][1] * emma.obj[[2]] +
          emma.obj[[1]][nrow(varcomp.values)] * diag(nrow(emma.obj[[2]]))

        ## Compute varcov matrix using var components from model
        vcovMatrix <- sommerFit$var.comp[1,1] * kinshipRed +
          diag(sommerFit$var.comp[nrow(sommerFit$var.comp), 1], nrow = nrow(kinshipRed))


      } else {
        ## if the kinship-matrix is the identity matrix ...
        vcovMatrix <- diag(sum(!is.na(gData$pheno[, trait])))
      }
    } else if (GLSMethod == 4) {
      vcovMatrix.list <- list()
      for (chr in gData$chromosomes) {
        chr.ind <- which(chr == gData$chromosomes)
        vcovMatrix.list[[chr.ind]] <- emma.obj.list[[chr.ind]][[1]][1] * emma.obj.list[[chr.ind]][[2]] +
          emma.obj.list[[chr.ind]][[1]][2] * diag(nrow(emma.obj.list[[chr.ind]][[2]]))
      }
    }



    # Name of the scan_GLS output-file
    output.file <- paste("results/",trait,".","output",suffix,".csv",sep="")

    ################################################
    # GWA USING scan_GLS (1), FaST-LMM (2), or R (3+4)
    maxScore <- max(gData$markers, na.rm = TRUE)

    ## Compute allele frequencies based on genotypes for which phenotypic data is available.
    allFreq <- colMeans(gData$markers[nonMissing, ])

    ## Compute variance of marker scores, based on genotypes for which phenotypic data is available.
    ## for inbreeders, this depends on maxScore. It is therefore scaled to marker scores 0, 1
    ## (or 0, 0.5, 1 if there are heterozygotes)
    allVar <- apply(gData$markers[nonMissing, ], 2, var) / maxScore ^ 2

    # allele frequencies based on all genotypes (trait-independent)
    allFreqTot <- colMeans(gData$markers)

    if (maxScore == 2) {
      allFreq <- allFreq / 2
      allFreqTot<- allFreqTot/ 2
    }

    if (!useMAF) {MAF <- MAC / length(nonMissing) - 0.00001}

    # if (GLSMethod ==1) {
    #
    #   GWAResult          <- scan_GLS(gData=gData,input.pheno=input.pheno,
    #     varcomp.file=varcomp.file,output.file=output.file,
    #     covariate.file=covariate.file,maf=MAF)
    #
    #   GWAResult          <- GWAResult[,setdiff(colnames(GWAResult),c('Nalleles','major_allele','minor_allele'))]
    #
    #   #######
    #   # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES
    #
    #   GC <- GenomicControlPvalues(pvals=GWAResult$pvalue,
    #     n.obs=sum(!is.na(gData$pheno[,trait])),
    #     n.cov=length(cov.cols))
    #
    #   inflationFactor  <- GC[[2]]
    #
    #   if (genomicControl) {
    #     GWAResult$pvalue <- GC[[1]]
    #   }
    #
    #   #allFreqSel
    #
    #   names(GWAResult)[which(names(GWAResult)=='stat')]   <- 'effect'
    #   names(GWAResult)[which(names(GWAResult)=='marker')] <- 'snp'
    #   names(GWAResult)[which(names(GWAResult)=='minorfreq')] <- 'allele.frequency'
    #
    #   # replace the minorfreq column by the frequencies calculated within R
    #   GWAResult$allele.frequency <- allFreq
    #
    #   # effects should be wrt the reference allele, not wrt the rare allele
    #   # (and The rare allele in scan_GLS is defined based on the whole panel; therefore
    #   #  the use of  allFreqTot)
    #   GWAResult$effect[allFreqTot > 0.5] <- -1 * GWAResult$effect[allFreqTot > 0.5]
    #
    #   # effects should be for a single allele, not for 2
    #   #if (maxScore==1) {GWAResult$effect <- 0.5 * GWAResult$effect}
    #   # REMOVED! already done in scan_GLS !
    #
    #   write.table(GWAResult,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
    # }
    #
    #
    #
    # if (GLSMethod ==2) {
    #
    #   if (kinship.type %in% 1:2) {
    #     gData$kinship <- gData$kinship / KinshipTransform(gData$kinship)
    #     WriteFastLmmKinship(gData=gData,trait=trait,file="output/temp_kinship.txt")
    #     FaST_LMM(gData=gData,trait.nr=trait,cov.cols=cov.cols,full.lm=full.scan,kinship.file="output/temp_kinship.txt",
    #       output.file=output.file,excludebyposition=excludebyposition,show.call=T,
    #       prefix=linux.prefix)
    #   }
    #
    #   markerMeans        <- apply(gData$markers[,nonMissing],1,mean)
    #   marker.sd           <- apply(gData$markers[,nonMissing],1,sd)
    #
    #   segMarkers <- which(markerMeans >= maxScore * MAF & (markerMeans <= maxScore * (1-MAF)))
    #
    #   GWAResult            <- ReadFastLmmResult(output.file)
    #   row.names(GWAResult) <- GWAResult$snp
    #
    #   if (allelic.substitution.effects) {
    #
    #     if (maxScore==1) {extra.constant <- 2} else {extra.constant <- 1}
    #
    #     GWAResult$SNPWeight   <- GWAResult$SNPWeight / (marker.sd[row.names(GWAResult)] * extra.constant)
    #     GWAResult$SNPWeightSE <- GWAResult$SNPWeightSE / (marker.sd[row.names(GWAResult)] * extra.constant)
    #
    #   }
    #
    #   GWAResult <- data.frame(GWAResult,
    #     R_LR_2=1 - exp(-2*(GWAResult$AltLogLike-GWAResult$NullLogLike)/GWAResult$N))
    #
    #
    #   if (nrow(GWAResult) < gData$N) {
    #
    #     non.segregating.marker.numbers <- setdiff(1:gData$N,segMarkers)
    #
    #     non.segregating.marker.names   <- setdiff(row.names(gData$markers),GWAResult$snp)
    #
    #     N.non.segregating              <- length(non.segregating.marker.names)
    #
    #     extra.block <- matrix(NA,N.non.segregating,ncol(GWAResult))
    #
    #     extra.block <- as.data.frame(extra.block)
    #
    #     colnames(extra.block) <- colnames(GWAResult)
    #
    #     extra.block[,c("chromosome","position")] <- gData$map[non.segregating.marker.names,c("chromosome","position")]
    #
    #     extra.block$snp <- as.character(non.segregating.marker.names)
    #
    #     rownames(extra.block) <- extra.block$snp
    #
    #     GWAResult <- rbind(GWAResult,extra.block)
    #
    #     GWAResult <- GWAResult[order(GWAResult$chromosome,GWAResult$position),]
    #
    #     GWAResult <- data.frame(GWAResult,allele.frequency=allFreq)
    #   }
    #
    #
    #   ########
    #   # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES
    #
    #   GC <- GenomicControlPvalues(pvals=GWAResult$pvalue,
    #     n.obs=sum(!is.na(gData$pheno[,trait])),
    #     n.cov=length(cov.cols))
    #
    #   inflationFactor  <- GC[[2]]
    #
    #   if (genomicControl) {
    #     GWAResult$pvalue <- GC[[1]]
    #   }
    #
    #   #######################################
    #
    #   colnames(GWAResult)[which(colnames(GWAResult)=='SNPWeight')]   <- 'effect'
    #   colnames(GWAResult)[which(colnames(GWAResult)=='SNPWeightSE')] <- 'effect_se'
    #
    #   ########################
    #
    #   GWAResult$pvalue[setdiff(1:gData$N,segMarkers)] <- NA
    #   GWAResult$effect[setdiff(1:gData$N,segMarkers)] <- NA
    #   GWAResult$effect_se[setdiff(1:gData$N,segMarkers)] <- NA
    #
    #   ##########################
    #
    #   #if (maxScore==1) {
    #   # effects should be wrt the reference allele, not wrt the rare allele
    #   # (and The rare allele in scan_GLS is defined based on the whole panel; therefore
    #   #  the use of  allFreqTot)
    #   GWAResult$effect[allFreqTot > 0.5] <- -1 * GWAResult$effect[allFreqTot > 0.5]
    #   #}
    #   #########################
    #
    #   write.table(GWAResult,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
    # }
    #
    #
    # #maxScore           <- max(gData$markers,na.rm=T)
    #
    #

    if (GLSMethod == 3) {
      Sigma <- vcovMatrix

      ## The following is based on the genotypes, not the replicates:
      markerMeans <- colMeans(gData$markers[nonMissing, ])
      segMarkers <- which(markerMeans >= maxScore * MAF & markerMeans <= maxScore * (1 - MAF))

      X.temp <- gData$markers[nonMissingRepId, segMarkers]
      rownames(X.temp) <- nonMissingRepId
      Y.temp <- phenoField[!is.na(phenoField[trait]), trait]
      names(Y.temp) <- nonMissingRepId

      GWAResult <- data.frame(snp = rownames(gData$map),
        gData$map,
        pValue = NA,
        effect = NA,
        effectSe = NA,
        RLR2 = NA,
        allFreq = allFreq,
        row.names = rownames(gData$map))

      if (sum(cov.cols) == 0) {
        GLSResult <- fastGLS(y = Y.temp, X = X.temp, Sigma = Sigma)
      } else {
        Z <- as.matrix(gData$pheno[!is.na(gData$pheno[, trait]), cov.cols])
        colnames(Z) <- names(gData$pheno)[cov.cols]
        rownames(Z) <- nonMissingRepId
        GLSResult <- fastGLSCov(y = Y.temp, X = X.temp, Sigma = Sigma, covs = Z)
      }
      GWAResult[segMarkers, c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult

      if (ncol(gData$genes) != 0) {
        GWAResult <- cbind(GWAResult, gene1 = gData$map$gene1, gene2 = gData$map$gene2)
      }

      ## Effects should be for a single allele, not for 2
      if (maxScore == 1) {
        GWAResult$effect <- 0.5 * GWAResult$effect
      }

      ## Calculate the genomic inflation factor and rescale p-values
      GC <- genomicControlPValues(pVals = GWAResult$pValue,
        nObs = sum(!is.na(phenoField[trait])),
        nCov = length(cov.cols))
      inflationFactor <- GC[[2]]
      if (genomicControl) {
        GWAResult$pvalue <- GC[[1]]
      }
      #write.table(GWAResult,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
      #rm(X.temp)
    } else if (GLSMethod ==4) {
      markerMeans <- rowMeans(gData$markers[, nonMissing])
      segMarkers <- which(markerMeans >= maxScore * MAF & markerMeans <= maxScore * (1 - MAF))

      Y.temp <- gData$pheno[!is.na(gData$pheno[, trait]), trait]
      colnames(Y.temp) <- nonMissingRepId

      GWAResult <- data.frame(snp = rownames(gData$map),
        gData$map,
        pValue = NA,
        effect = NA,
        effectSe = NA,
        RLR2 = NA,
        allFreq = allFreq,
        row.names = rownames(gData$map))

      if (sum(cov.cols) == 0) {
        for (chr in gData$chromosomes) {
          X.temp <- t(gData$markers[intersect(which(gData$map$chromosome == chr), segMarkers),
            nonMissingRepId])
          rownames(X.temp) <- nonMissingRepId
          Sigma <- vcovMatrix.list[[which(chr == gData$chromosomes)]]
          GLSResult <- fastGLS(y = Y.temp, X = X.temp, Sigma = Sigma)
        }
      } else {
        Z <- as.matrix(gData$pheno[!is.na(gData$pheno[,trait]), cov.cols])
        colnames(Z) <- names(gData$pheno)[cov.cols]
        rownames(Z)<- nonMissingRepId
        for (chr in gData$chromosomes) {
          X.temp <- t(gData$markers[intersect(which(gData$map$chromosome == chr), segMarkers),
            nonMissingRepId])
          row.names(X.temp) <- nonMissingRepId
          Sigma <- vcovMatrix.list[[which(chr == gData$chromosomes)]]
          GLSResult <- fastGLSCov(y = Y.temp, X = X.temp, Sigma = Sigma, covs = Z)
        }
      }
      GWAResult[intersect(which(gData$map$chromosome == chr), segMarkers),
        c("pValue", "effect", "effectSe", "RLR2")] <- GLSResult

      ## Effects should be for a single allele, not for 2
      if (maxScore==1) {
        GWAResult$effect <- 0.5 * GWAResult$effect
      }

      ## Calculate the genomic inflation factor and rescale p-values
      GC <- genomicControlPValues(pVals = GWAResult$pValue,
        nObs = sum(!is.na(phenoField[trait])),
        nCov = length(cov.cols))
      inflationFactor <- GC[[2]]
      if (genomicControl) {
        GWAResult$pValue <- GC[[1]]
      }
      #write.table(GWAResult,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
      #rm(X.temp)
    }

    nEff <- sum(!is.na(GWAResult$pValue))
    ## Calculate the significance threshold
    ## When boundType is 1, 3 or 4, determine the LOD-threshold

    if (boundType == 1) {
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
        snpStatus <- rep("significant snp", length(signSnp))
      }
      if (GLSMethod==1) {
        #   effects                  <- GWAResult$effect[snpSelection]
        #
        #   snpVar     <- 4 * effects^2 * allVar[snpSelection]
        #   propSnpVar <- snpVar / var(as.numeric(na.omit(gData$pheno[,trait])))
        #
        #   snp.sig             <- cbind(snp.sig,data.frame(reference.allele.freq=allFreqSel,
        #     effect.size=effects,
        #     propSnpVar=propSnpVar))
        #   #perc.of.phenetic.var=explained.phenotypic.variance))
        # } else if (GLSMethod==2) {
        #   effects                  <- GWAResult$effect[snpSelection]
        #   #snpVar     <- 4 * effects^2 * allFreqSel * (1 - allFreqSel)
        #   #explained.phenotypic.variance <- 100 * snpVar / var(as.numeric(na.omit(gData$pheno[,trait])))
        #   #
        #   snp.sig             <- cbind(snp.sig,data.frame(reference.allele.freq=allFreqSel,
        #     effect.size=effects,
        #     proportion.var.LRT=GWAResult$R_LR_2[snpSelection]))
        #   #                             perc.of.phenetic.var=explained.phenotypic.variance))
      } else if (GLSMethod %in% 3:4) {
        effects <- GWAResult$effect[snpSelection]
        snpVar <- 4 * effects ^ 2 * allVar[snpSelection]
        propSnpVar <- snpVar / as.numeric(var(na.omit(phenoField[trait])))
      }

      snpSig <- data.frame(marker = colnames(gData$markers)[snpSelection],
        gData$map[snpSelection, ],
        pValue = GWAResult$pValue[snpSelection],
        snpStatus,
        allFreq = allFreq[snpSelection],
        effects,
        propVarLRT = GWAResult$RLR2[snpSelection],
        propSnpVar = propSnpVar,
        stringsAsFactors = FALSE)


      #############

      # restore for arabidopsis ?
      #accessions.with.the.rare.allele <- rep("",nrow(snp.sig))
      #for (snp in 1:nrow(snp.sig)) {
      #  if (snp.sig$reference.allele.freq[snp] < 0.5) {rare.allele<-1} else {rare.allele<-0}
      #  temp.obj <- gData$markers[snp.sig$marker[snp],unique(gData$pheno$genotype[!is.na(gData$pheno[,trait])])]
      #  temp.names <- names(temp.obj)
      #  accessions.with.the.rare.allele[snp] <- paste(temp.names[temp.obj==rare.allele],collapse=",")
      #}
      #snp.sig <- cbind(snp.sig,data.frame(accessions.with.the.rare.allele=accessions.with.the.rare.allele))

      # creating this file
      # snp.output.file       <- paste("results/","significant.snps.",trait,suffix,".csv",sep="")

      #write.table(snp.sig,file=snp.output.file,quote=FALSE,row.names=FALSE,sep=",")

      #if (make.all.trait.significant.snp.file) {
      #  snp.sig <- data.frame(snp.sig[,1:(ncol(snp.sig)-1)],trait.name=rep(trait,nrow(snp.sig)),snp.sig[,ncol(snp.sig)])
      #  write.table(snp.sig,file=paste("results/",all.trait.significant.snp.file,sep=""),append=T,quote=FALSE,row.names=FALSE,col.names=F,sep=",")
      #}


    }

    # if (make.all.trait.h2.file) {
    #   all.traits$h2[which(trait.numbers==trait)] <- (varcomp.values[1,1])/(sum(varcomp.values[,1])) # which(names(gData$pheno)==trait)-min(trait.numbers)+1
    #   all.traits$inflation[which(trait.numbers==trait)] <- inflationFactor
    #   if (kinship.type %in% 3:4) {all.traits$n.markers.in.K[which(trait.numbers==trait)] <- testK$m.optimal}
    # }
    #
    # if (make.one.file.with.all.lod.scores) {
    #   all.lod.scores.frame[,which(names(all.lod.scores.frame)==trait)] <- -log10(ReplaceNaByOne(GWAResult$pvalue))
    # }
    #
    # if (make.one.file.with.all.snp.effects) {
    #   if (GLSMethod==1) {all.snp.effects.frame[,which(names(all.snp.effects.frame)==trait)] <- GWAResult$effect}
    #   if (GLSMethod %in% 2:4) {all.snp.effects.frame[,which(names(all.snp.effects.frame)==trait)] <- GWAResult$effect}
    # }

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


    ##################################################################################################

    # Make a plot of the LOD-profile: (pdf and jpeg)
    # LODThr.plot  <- 0
    # if (boundType %in% 1:4) {LODThr.plot <- LODThr}
    # if (boundType == 5) {if (!no.significant.snps) {LODThr.plot <- min(-log10(GWAResult$pvalue[snpSelection]))}}
    # #
    # x.effects=integer(0)
    # effects.size=numeric(0)
    #
    # if ("real.effects" %in% names(gData)) {
    #   if (nrow(gData$real.effects$locations)>0) {
    #     x.effects=gData$real.effects$locations[,trait-1]
    #     effects.size=gData$real.effects$sizes[,trait-1]
    #   }
    # }
    #
    # if (no.significant.snps) {
    #   x.sig <- integer(0)
    # } else {
    #   x.sig <- which(snpSelection)[which(snp.sig$snpStatus=="significant snp")]
    # }
    #
    #
    # MakeLodPlotWithChromosomeColors(xvalues=gData$map$cum.position,yvalues=-log10(ReplaceNaByOne(GWAResult$pvalue)),plot.title=trait,gData=gData,
    #   col.palette = rep(c("royalblue","maroon"),50)[1:gData$nchr],file.name=paste("plots/",trait,suffix,".jpeg",sep=""),
    #   x.sig=x.sig,chr.boundaries=cumsum(gData$chr.lengths.bp)[-1],y.thr=LODThr.plot,
    #   x.effects=x.effects,effects.size=effects.size,plot.type="d")
    #
    # if (!jpeg.only) {
    #   MakeLodPlotWithChromosomeColors(xvalues=gData$map$cum.position,yvalues=-log10(GWAResult$pvalue),plot.title=trait,gData=gData,
    #     col.palette = rep(c("royalblue","maroon"),50)[1:gData$nchr],
    #     file.name=paste("plots/",trait,suffix,".pdf",sep=""),x.sig=x.sig,
    #     chr.boundaries=cumsum(gData$chr.lengths.bp)[-1],y.thr=LODThr.plot,x.effects=x.effects,
    #     effects.size=effects.size,jpeg.plot=FALSE,plot.type="d")
    # }
    #
    # QQplotPvalues2(pvalues=GWAResult$pvalue,main.title=trait,file.name=paste0("plots/QQplot_",trait,suffix,".png"))
    #

  }   # end for (trait in trait.numbers)


  ########################################################################################################################
  ########################################################################################################################

  # if (make.all.trait.h2.file) {
  #   write.csv(all.traits,file=paste("results/",all.trait.h2.file,sep=""),row.names=F)
  # }
  #
  # if (make.one.file.with.all.lod.scores) {
  #   write.csv(all.lod.scores.frame,file=paste("results/",file.with.all.lod.scores,sep=""),quote=F)
  # }
  #
  # if (make.one.file.with.all.snp.effects) {
  #   write.csv(all.snp.effects.frame,file=paste("results/",file.with.all.snp.effects,sep=""),quote=F)
  # }

  #if (covariables) {
  #  # After the loop over traits
  #  gData$pheno <- temp.pheno
  #}
  return(GWAResult = GWAResult, signSnp = signSnp)
}
