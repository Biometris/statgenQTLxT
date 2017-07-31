#' run single trait GWAS

#' @param gData an object of class \code{gData}
#' @param GLSMethod an integer indicating the software used to estimate the marker effects.
#' \enumerate{
#' \item{scan_GLS (not available when there are heterozygotes) - NOT IMPLEMENTED YET}
#' \item{Fast-LMM - NOT IMPLEMENTED YET}
#' \item{within R}
#' \item{within R, with chromosome specific kinship matrices; similar to the approach of Rincent et al.
#' (genetics, 2014)}
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
#'
#' @references Rincent et al. (2014) ....






runSingleTraitGwas <- function (gData,
  GLSMethod = 3,
  sizeIncludedRegion = 0,
  minR2,
  useMAF = TRUE,
  MAF = 0.05,
  MAC = 10) {


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




  if (strsplit(data.path,split=" ")[[1]][1]!=data.path) {
    stop("ERROR: data.path should not contain spaces.")
  }
  if (strsplit(script.path,split=" ")[[1]][1]!=script.path) {
    stop("ERROR: scipt.path should not contain spaces.")
  }

  if (useMAF) {
    if (MAF == 0) {MAF <- 1e-6}
  } else {
    if (MAC == 0) {MAC <- 1}
  }

  setwd(data.path)

  # Create the subdirectory structure. If it already exists, nothing is done.
  dir.create("results", showWarnings = F)
  dir.create("plots", showWarnings = F)
  dir.create("output", showWarnings = F)

  ####################################################################################################################################################
  # load functions and data:

  if (GLSMethod==4) {
    if (!(all(paste0('kchr',gData$chromosomes) %in% names(gData)))) {
      stop("If GLSMethod=4, gData should contain chromosome-specific-kinship matrices kchr1, kchr2, etc.")
    }
    if (kinship.type==2) {
      kinship.type <- 1
      cat(paste0('Warning: GLSMethod=4: Therefore the kinship matrix contained in ',alternative.kinship.name,' is ignored.','\n'))
    }
    if (reml.algo=="asreml") {
      reml.algo <- 'emma'
      cat(paste0('Warning: GLSMethod=4: Therefore reml.algo is set to emma','\n'))
    }
    if (h2.fixed==TRUE) {
      h2.fixed <- FALSE
      cat(paste0('Warning: GLSMethod=4: Therefore h2.fixed is set to FALSE','\n'))
    }
    #if (!all(gData$chromosomes==1:gData$nchr)) {
    #  stop("When GLSMethod=4, the chromosomes in gData (indicated by gData$chromosomes)\n should be labeled 1,...,k, where k is the number given by gData$nchr")
    #}
    if (make.all.trait.h2.file==TRUE) {
      make.all.trait.h2.file <- FALSE
      cat(paste0('Warning: GLSMethod=4: Therefore make.all.trait.h2.file is set to FALSE','\n'))
    }
  }

  source(paste0(script.path,"gwas_functions.R"))

  if (add.phenotypic.data) {
    gData <- AddPhenoData(gData=gData,csv.file.name=csv.file.name,add.var.means=F,make.pheno.image=F,pheno.image.name="",
      which.columns.as.factor=which.columns.as.factor)

    if (length(which.columns.as.factor) > 0) {
      if (!all(which.columns.as.factor %in% 2:ncol(gData$pheno))) {
        stop("ERROR: which.columns.as.factor should be integer(0) (the default), or a number \n between 2 and the number of columns in the .csv file with the phenotypes")
      }
    }

    for (k in which.columns.as.factor) {
      if ("NA" %in% levels(gData$pheno[,k])) {
        #gData$pheno[which(is.na(gData$pheno[,2])),2] <- "NA"
        is.na(gData$pheno[,k])[gData$pheno[,k]=="NA"] <- TRUE
        suppressWarnings(gData$pheno[,k] <- factor(gData$pheno[,k],exclude="NA"))
      }
    }
  }

  if (!all(trait.numbers %in% setdiff(2:ncol(gData$pheno),which.columns.as.factor))) {
    stop("ERROR: all numbers in the vector trait.numbers should be between 2 and the \n number of columns in the .csv file with the phenotypes. \n In addition, colmun-numbers declared in which.columns.as.factor cannot be in trait.numbers")
  }

  ##################################################################################################################

  # install and load R-packages asreml and multtest (asreml needs to be installed from a zip file)
  if (!(reml.algo %in% c("asreml","emma"))) {reml.algo <- "emma"}
  if (kinship.type==4) {reml.algo     <- "emma"}

  if (reml.algo %in% c("asreml")) {
    if (!is.installed("asreml")) {cat("ERROR: first install asreml, or use emma","\n")}
    library(asreml)
  } #else {
  #  source(paste(script.path,"emma.R",sep=""))
  #}

  if (var.explained.snp.random & (!is.installed("asreml"))) {
    cat("WARNING: asreml is not available; the option var.explained.snp.random is set to FALSE","\n")
    var.explained.snp.random <- FALSE
  }


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
    #genomic.control <- T
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


  if (GLSMethod==2) {
    if (!file.exists(paste(gData$external$plink.name,".bed",sep=""))) {
      cat("Creating PLINK binary files; this may take a while.\n")
      MakePlinkFiles(gData,file.name=plink.file.name)
      gData <- AddBinaryPlinkFilesToGwasObj(gData=gData,plink.name=plink.file.name)
      save(gData,file=r.image.name)
      cat("",file="fastlmmLog.txt")
      #   if (show.call) {cat(command.string,"\n",file="fastlmmLog.txt")}
    }
    #!# minor allele fr's should be based on only those accessions for which phenotypic values are available
    #minor.allele.frequencies <- apply(gData$markers,1,function(x){mean(as.numeric(x))})
    #minor.allele.frequencies[minor.allele.frequencies>.5] <-  1-minor.allele.frequencies[minor.allele.frequencies>.5]

  }

  if (reml.algo == "asreml") {
    gData$kinship.asreml <- MakeKinshipAsreml(gData$kinship,genotype.names=gData$plant.names)
  }

  ################################################################
  ####### covariates

  # if no covariates are to be used, cov.frame will be an dataframe with one column ("mu": vector of ones).
  # Otherwise, it will be a dataframe with ones in the first
  # column, and the actual covariates in subsequent columns

  if (!covariables) {
    cov.cols  <- integer(0)
    covariate.file  <- ""
  }

  if (covariables) {
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

  if (make.all.trait.significant.snp.file) {cat("",file=paste("results/",all.trait.significant.snp.file,sep=""))}

  if (make.all.trait.h2.file) {
    all.traits <- data.frame(trait=names(gData$pheno)[trait.numbers],h2=rep(NA,length(trait.numbers)),inflation=rep(NA,length(trait.numbers)))
  }

  if (make.one.file.with.all.lod.scores) {
    all.lod.scores.frame <- data.frame(gData$map,matrix(NA,gData$N,length(trait.numbers)))
    names(all.lod.scores.frame)[-(1:ncol(gData$map))] <- names(gData$pheno)[trait.numbers]
  }

  if (make.one.file.with.all.snp.effects) {
    all.snp.effects.frame <- data.frame(gData$map,matrix(NA,gData$N,length(trait.numbers)))
    names(all.snp.effects.frame)[-(1:ncol(gData$map))] <- names(gData$pheno)[trait.numbers]
  }

  ###################################################################################################


  if (covariables) {
    temp.pheno <- gData$pheno
  }

  #######################################################################################################################################################################
  #######################################################################################################################################################################
  #######################################################################################################################################################################

  for (tr.n in trait.numbers) {  # loop over traits
    #tr.n = trait.numbers[1]

    start.date <- date()

    if (strsplit(names(gData$pheno)[tr.n],split=" ")[[1]][1]!=names(gData$pheno)[tr.n]) {
      stop(paste("ERROR: the variable name",names(gData$pheno)[tr.n],"should not contain spaces."))
    }

    if (covariables) {
      # if covariates are used, all phenotypic values (for the current trait) for which at least one covariate is missing, are set to missing
      gData$pheno[apply(as.data.frame(gData$pheno[,cov.cols]),1,function(x){sum(is.na(x))})>0,tr.n]  <- NA
    }


    # the name of the trait, which will be the basis of many file names:
    trait       <- names(gData$pheno)[tr.n]


    ###############################################################
    # ESTIMATION OF THE VARIANCE COMPONENTS WITH EITHER EMMA OR ASREML


    if (reml.algo=="asreml") {

      random.formula <- as.formula("~ giv(genotype,var=T)")

      if (sum(cov.cols)!=0) {  # reml.formula is the formula for the fixed part of the model
        reml.formula        <- as.formula(paste(paste(names(gData$pheno)[tr.n],"~"),paste(names(gData$pheno)[cov.cols],collapse="+")))
      } else {
        reml.formula        <- as.formula(paste(names(gData$pheno)[tr.n],"~ 1"))
      }

      reml.obj            <- asreml(aom=F,maxiter = 25,fixed= reml.formula,data=gData$pheno[!is.na(gData$pheno[,tr.n]),],
        random = random.formula,na.method.X="omit",ginverse = list(genotype=gData$kinship.asreml))

      varcomp.values      <- data.frame(var.comp.values=summary(reml.obj)$varcomp$component)
      varcomp.std         <- summary(reml.obj)$varcomp$std.error

    }

    if (reml.algo=="emma" & GLSMethod!=4) {
      if (h2.fixed) {

        kinship.reduced  <- MakeScanGlsKinship(1,0,gData$kinship,plant.names=gData$plant.names,gData$pheno,tr.n=tr.n)
        Knr <- kinship.reduced / KinshipTransform(kinship.reduced)
        Knr <- kinship.reduced

        delta <- h2 / (1-h2)
        evk <- eigen(Knr)
        U <- evk[[2]]
        v <- evk[[1]]
        xbar <- apply(U,2,sum)
        ybar <- as.numeric(t(U) %*% as.matrix(gData$pheno[!is.na(gData$pheno[,tr.n]),tr.n]))
        sbar <- 1/(v+delta)
        betahat <- sum(sbar * xbar * ybar) / sum(sbar * xbar * xbar)
        sigmaG2 <- mean(sbar * (ybar - betahat * xbar)^2)
        sigmaE2 <- delta * sigmaG2
        varcomp.values      <- data.frame(var.comp.values=c(sigmaG2/KinshipTransform(kinship.reduced),sigmaE2))

      } else {

        emma.obj            <- RunEmma(gData=gData,tr.n,cov.cols)
        varcomp.values      <- data.frame(var.comp.values=emma.obj[[1]])

      }

      varcomp.std         <- c(NA,NA)

    }

    #emma.obj = RunEmma(gData=gData,tr.n,cov.cols)
    #vr = emma.obj[[1]][1] * emma.obj[[2]] + emma.obj[[1]][2] * diag(nrow(emma.obj[[2]]))
    # sum(abs(vcov.matrix - vr))

    #!# to do : if the kinship-matrix is the identity matrix ...


    if (reml.algo=="emma" & GLSMethod==4) {

      emma.obj.list       <- list()
      varcomp.values.list <- list()

      for (chr in gData$chromosomes) {

        emma.obj.list[[which(chr==gData$chromosomes)]]  <- RunEmma(gData=gData,tr.n,cov.cols,
          K.user=gData[[which(names(gData)==paste0('kchr',chr))]])
        varcomp.values.list[[which(chr==gData$chromosomes)]]      <- data.frame(var.comp.values=emma.obj.list[[which(chr==gData$chromosomes)]][[1]])

      }

      #varcomp.std         <- c(NA,NA)

    }


    ##########################################################
    # CREATE PHENO-TYPE FILE + VARCOMP-FILE FOR  SCAN-GLS / FAST-LMM
    #

    if (GLSMethod==1) {

      varcomp.file        <- paste("output/",trait,".","varcomp",".csv",sep="")

      # Name of the phenotype-file for scan_GLS:
      input.pheno <- paste("output/",trait,".csv",sep="")

      if (sum(abs(gData$kinship-diag(nrow(gData$kinship))))>0.0001 ) { # if the kinship-matrix is (numerically) different from the identity matrix ...
        # to do : the case where the kinship matrix is diagonal, but has unequal variances
        vcov.matrix <- MakeScanGlsKinship(varcomp.values[1,],varcomp.values[nrow(varcomp.values),],
          gData$kinship,plant.names=gData$plant.names,
          gData$pheno,tr.n=tr.n)
      } else {  # if the kinship-matrix is the identity matrix ...
        vcov.matrix <- diag(sum(!is.na(gData$pheno[,tr.n])))
      }
      write.table(vcov.matrix,file=varcomp.file,sep=",",quote=F,col.names=F,row.names=F)
      MakePhenoFile(pheno.object=gData$pheno,col.number=tr.n,file.name=input.pheno)
    }

    if (GLSMethod==3) {
      if (sum(abs(gData$kinship-diag(nrow(gData$kinship))))>0.0001 ) { # if the kinship-matrix is (numerically) different from the identity matrix ...
        # to do : the case where the kinship matrix is diagonal, but has unequal variances
        vcov.matrix <- emma.obj[[1]][1] * emma.obj[[2]] + emma.obj[[1]][nrow(varcomp.values)] * diag(nrow(emma.obj[[2]]))
      } else {  # if the kinship-matrix is the identity matrix ...
        vcov.matrix <- diag(sum(!is.na(gData$pheno[,tr.n])))
      }
    }

    if (GLSMethod==4) {
      vcov.matrix.list <- list()
      for (chr in gData$chromosomes) {
        chr.ind <- which(chr==gData$chromosomes)
        vcov.matrix.list[[chr.ind]] <- emma.obj.list[[chr.ind]][[1]][1] * emma.obj.list[[chr.ind]][[2]] + emma.obj.list[[chr.ind]][[1]][2] * diag(nrow(emma.obj.list[[chr.ind]][[2]]))
      }
    }

    ################################################

    non.missing     <- unique(gData$pheno$genotype[!is.na(gData$pheno[,tr.n])])
    kinship.reduced <- gData$kinship[non.missing,non.missing]

    # new:
    non.missing.rep.id <- row.names(gData$pheno[!is.na(gData$pheno[,tr.n]),])
    non.missing.rep.geno <- gData$pheno$genotype[!is.na(gData$pheno[,tr.n])]

    # Name of the scan_GLS output-file
    output.file         <- paste("results/",trait,".","output",suffix,".csv",sep="")

    ################################################
    # GWA USING scan_GLS (1), FaST-LMM (2), or R (3+4)

    max.score           <- max(gData$markers,na.rm=T)

    # allele frequencies based on genotypes for which phenotypic data is available (trait-dependent)
    all.allele.frequencies    <- apply(gData$markers[,unique(gData$pheno$genotype[!is.na(gData$pheno[,tr.n])])],1,mean)

    # variances of marker scores, based on genotypes for which phenotypic data is available (trait-dependent)
    # for inbreeders, this depends on max.score. It is therefore scaled to marker scores 0,1 (or 0,0.5,1 if there are heterozygotes)
    all.allele.variances      <- apply(gData$markers[,unique(gData$pheno$genotype[!is.na(gData$pheno[,tr.n])])],1,var) / (max.score)^2

    # allele frequencies based on all genotypes (trait-independent)
    all.allele.frequencies2   <- apply(gData$markers,1,mean)

    if (max.score==2) {all.allele.frequencies <- all.allele.frequencies / 2}

    if (max.score==2) {all.allele.frequencies2<- all.allele.frequencies2/ 2}

    if (useMAF==FALSE) {MAF <- MAC / length(non.missing) - 0.00001}

    if (GLSMethod ==1) {

      GWA.result          <- scan_GLS(gData=gData,input.pheno=input.pheno,
        varcomp.file=varcomp.file,output.file=output.file,
        covariate.file=covariate.file,maf=MAF)

      GWA.result          <- GWA.result[,setdiff(colnames(GWA.result),c('Nalleles','major_allele','minor_allele'))]

      #######
      # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES

      GC <- GenomicControlPvalues(pvals=GWA.result$pvalue,
        n.obs=sum(!is.na(gData$pheno[,tr.n])),
        n.cov=length(cov.cols))

      inflation.factor  <- GC[[2]]

      if (genomic.control) {
        GWA.result$pvalue <- GC[[1]]
      }

      #allele.frequencies

      names(GWA.result)[which(names(GWA.result)=='stat')]   <- 'effect'
      names(GWA.result)[which(names(GWA.result)=='marker')] <- 'snp'
      names(GWA.result)[which(names(GWA.result)=='minorfreq')] <- 'allele.frequency'

      # replace the minorfreq column by the frequencies calculated within R
      GWA.result$allele.frequency <- all.allele.frequencies

      # effects should be wrt the reference allele, not wrt the rare allele
      # (and The rare allele in scan_GLS is defined based on the whole panel; therefore
      #  the use of  all.allele.frequencies2)
      GWA.result$effect[all.allele.frequencies2 > 0.5] <- -1 * GWA.result$effect[all.allele.frequencies2 > 0.5]

      # effects should be for a single allele, not for 2
      #if (max.score==1) {GWA.result$effect <- 0.5 * GWA.result$effect}
      # REMOVED! already done in scan_GLS !

      write.table(GWA.result,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
    }



    if (GLSMethod ==2) {

      if (kinship.type %in% 1:2) {
        gData$kinship <- gData$kinship / KinshipTransform(gData$kinship)
        WriteFastLmmKinship(gData=gData,tr.n=tr.n,file="output/temp_kinship.txt")
        FaST_LMM(gData=gData,trait.nr=tr.n,cov.cols=cov.cols,full.lm=full.scan,kinship.file="output/temp_kinship.txt",
          output.file=output.file,excludebyposition=excludebyposition,show.call=T,
          prefix=linux.prefix)
      }

      marker.means        <- apply(gData$markers[,non.missing],1,mean)
      marker.sd           <- apply(gData$markers[,non.missing],1,sd)

      segregating.markers <- which(marker.means >= max.score * MAF & (marker.means <= max.score * (1-MAF)))

      GWA.result            <- ReadFastLmmResult(output.file)
      row.names(GWA.result) <- GWA.result$snp

      if (allelic.substitution.effects) {

        if (max.score==1) {extra.constant <- 2} else {extra.constant <- 1}

        GWA.result$SNPWeight   <- GWA.result$SNPWeight / (marker.sd[row.names(GWA.result)] * extra.constant)
        GWA.result$SNPWeightSE <- GWA.result$SNPWeightSE / (marker.sd[row.names(GWA.result)] * extra.constant)

      }

      GWA.result <- data.frame(GWA.result,
        R_LR_2=1 - exp(-2*(GWA.result$AltLogLike-GWA.result$NullLogLike)/GWA.result$N))


      if (nrow(GWA.result) < gData$N) {

        non.segregating.marker.numbers <- setdiff(1:gData$N,segregating.markers)

        non.segregating.marker.names   <- setdiff(row.names(gData$markers),GWA.result$snp)

        N.non.segregating              <- length(non.segregating.marker.names)

        extra.block <- matrix(NA,N.non.segregating,ncol(GWA.result))

        extra.block <- as.data.frame(extra.block)

        colnames(extra.block) <- colnames(GWA.result)

        extra.block[,c("chromosome","position")] <- gData$map[non.segregating.marker.names,c("chromosome","position")]

        extra.block$snp <- as.character(non.segregating.marker.names)

        rownames(extra.block) <- extra.block$snp

        GWA.result <- rbind(GWA.result,extra.block)

        GWA.result <- GWA.result[order(GWA.result$chromosome,GWA.result$position),]

        GWA.result <- data.frame(GWA.result,allele.frequency=all.allele.frequencies)
      }


      ########
      # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES

      GC <- GenomicControlPvalues(pvals=GWA.result$pvalue,
        n.obs=sum(!is.na(gData$pheno[,tr.n])),
        n.cov=length(cov.cols))

      inflation.factor  <- GC[[2]]

      if (genomic.control) {
        GWA.result$pvalue <- GC[[1]]
      }

      #######################################

      colnames(GWA.result)[which(colnames(GWA.result)=='SNPWeight')]   <- 'effect'
      colnames(GWA.result)[which(colnames(GWA.result)=='SNPWeightSE')] <- 'effect_se'

      ########################

      GWA.result$pvalue[setdiff(1:gData$N,segregating.markers)] <- NA
      GWA.result$effect[setdiff(1:gData$N,segregating.markers)] <- NA
      GWA.result$effect_se[setdiff(1:gData$N,segregating.markers)] <- NA

      ##########################

      #if (max.score==1) {
      # effects should be wrt the reference allele, not wrt the rare allele
      # (and The rare allele in scan_GLS is defined based on the whole panel; therefore
      #  the use of  all.allele.frequencies2)
      GWA.result$effect[all.allele.frequencies2 > 0.5] <- -1 * GWA.result$effect[all.allele.frequencies2 > 0.5]
      #}
      #########################

      write.table(GWA.result,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
    }


    #max.score           <- max(gData$markers,na.rm=T)



    if (GLSMethod ==3) {
      Sigma <- vcov.matrix
      #rownames(Sigma) <-  non.missing.rep.id
      #colnames(Sigma) <-  non.missing.rep.id

      # the following is based on the genotypes, not the replicates:

      marker.means        <- apply(gData$markers[,non.missing],1,mean)
      segregating.markers <- which(marker.means >= max.score * MAF & (marker.means <= max.score * (1-MAF)))

      X.temp            <- t(gData$markers[segregating.markers,non.missing.rep.geno])
      row.names(X.temp) <- non.missing.rep.id

      Y.temp              <- gData$pheno[!is.na(gData$pheno[,tr.n]),tr.n]
      names(Y.temp)       <- non.missing.rep.id

      GWA.result          <- data.frame(gData$map[,1:2],pvalue=rep(NA,gData$N),
        effect=rep(NA,gData$N),
        effect_se=rep(NA,gData$N),
        R_LR_2=rep(NA,gData$N))

      row.names(GWA.result) <-   row.names(gData$markers)

      if (sum(cov.cols)==0) {
        ww <- fastGLS_with_effect_sizes(Y=Y.temp,X=X.temp,Sigma=Sigma)
        GWA.result$pvalue[segregating.markers]    <- ww$pvalue
        GWA.result$effect[segregating.markers]    <- ww$beta
        GWA.result$effect_se[segregating.markers] <- ww$beta_se
        GWA.result$R_LR_2[segregating.markers]    <- ww$R_LR_2
      } else {
        Z <- as.matrix(gData$pheno[!is.na(gData$pheno[,tr.n]),cov.cols])
        colnames(Z) <- names(gData$pheno)[cov.cols]
        rownames(Z)<- non.missing.rep.id
        GWA.result$pvalue[segregating.markers]          <- fastGLS_cof(Y=Y.temp,X=X.temp,cofs=Z,Sigma=Sigma)

        ww <- fastGLS_cof_with_effect_sizes(Y=Y.temp,X=X.temp,cofs=Z,Sigma=Sigma)
        GWA.result$pvalue[segregating.markers] <- ww$pvalue
        GWA.result$effect[segregating.markers] <- ww$beta
        GWA.result$effect_se[segregating.markers] <- ww$beta_se
        GWA.result$R_LR_2[segregating.markers] <- ww$R_LR_2
      }

      GWA.result <- data.frame(GWA.result,allele.frequency=all.allele.frequencies)

      if (ncol(gData$genes)!=0) {
        old.names           <- colnames(GWA.result)
        GWA.result          <- cbind(GWA.result,gene1=gData$map$gene1,gene2=gData$map$gene2)
        colnames(GWA.result)<- c(old.names,"gene1","gene2")
      }

      # effects should be for a single allele, not for 2
      if (max.score==1) {GWA.result$effect <- 0.5 * GWA.result$effect}

      ########
      # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES

      GC <- GenomicControlPvalues(pvals=GWA.result$pvalue,
        n.obs=sum(!is.na(gData$pheno[,tr.n])),
        n.cov=length(cov.cols))

      inflation.factor  <- GC[[2]]

      if (genomic.control) {
        GWA.result$pvalue <- GC[[1]]
      }
      #######################################
      # already done
      #if (max.score==1) {
      #  GWA.result$effect   <- 0.5 * GWA.result$effect
      #}

      #############################

      GWA.result <- data.frame(snp=rownames(gData$markers),GWA.result)


      ########################
      write.table(GWA.result,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
      rm(X.temp)
    }



    if (GLSMethod ==4) {

      marker.means        <- apply(gData$markers[,non.missing],1,mean)
      segregating.markers <- which(marker.means >= max.score * MAF & (marker.means <= max.score * (1-MAF)))

      Y.temp              <- gData$pheno[!is.na(gData$pheno[,tr.n]),tr.n]
      names(Y.temp)       <- non.missing.rep.id

      GWA.result          <- data.frame(gData$map[,1:2],pvalue=rep(NA,gData$N),
        effect=rep(NA,gData$N),
        effect_se=rep(NA,gData$N),
        R_LR_2=rep(NA,gData$N))

      row.names(GWA.result) <-   row.names(gData$markers)

      if (sum(cov.cols)==0) {

        for (chr in gData$chromosomes) {

          X.temp            <- t(gData$markers[intersect(which(gData$map$chromosome==chr),segregating.markers),non.missing.rep.geno])
          row.names(X.temp) <- non.missing.rep.id
          Sigma <- vcov.matrix.list[[which(chr==gData$chromosomes)]]

          ww <- fastGLS_with_effect_sizes(Y=Y.temp,X=X.temp,Sigma=Sigma)
          GWA.result$pvalue[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$pvalue
          GWA.result$effect[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$beta
          GWA.result$effect_se[intersect(which(gData$map$chromosome==chr),segregating.markers)] <- ww$beta_se
          GWA.result$R_LR_2[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$R_LR_2
        }
      } else {

        Z <- as.matrix(gData$pheno[!is.na(gData$pheno[,tr.n]),cov.cols])
        colnames(Z) <- names(gData$pheno)[cov.cols]
        rownames(Z)<- non.missing.rep.id

        for (chr in gData$chromosomes) {

          X.temp            <- t(gData$markers[intersect(which(gData$map$chromosome==chr),segregating.markers),non.missing.rep.geno])
          row.names(X.temp) <- non.missing.rep.id
          Sigma <- vcov.matrix.list[[which(chr==gData$chromosomes)]]
          #GWA.result$pvalue[segregating.markers]          <- fastGLS_cof(Y=Y.temp,X=X.temp,cofs=Z,Sigma=Sigma)

          ww <- fastGLS_cof_with_effect_sizes(Y=Y.temp,X=X.temp,cofs=Z,Sigma=Sigma)
          GWA.result$pvalue[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$pvalue
          GWA.result$effect[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$beta
          GWA.result$effect_se[intersect(which(gData$map$chromosome==chr),segregating.markers)] <- ww$beta_se
          GWA.result$R_LR_2[intersect(which(gData$map$chromosome==chr),segregating.markers)]    <- ww$R_LR_2
        }
      }

      # effects should be for a single allele, not for 2
      if (max.score==1) {GWA.result$effect <- 0.5 * GWA.result$effect}

      ######################################################################################
      # CALCULATE THE GENOMIC INFLATION FACTOR, AND, IF THIS OPTION IS CHOSEN, RESCALE THE P-VALUES

      GC <- GenomicControlPvalues(pvals=GWA.result$pvalue,
        n.obs=sum(!is.na(gData$pheno[,tr.n])),
        n.cov=length(cov.cols))

      inflation.factor  <- GC[[2]]

      if (genomic.control) {
        GWA.result$pvalue <- GC[[1]]
      }

      #######################################

      GWA.result <- data.frame(snp=rownames(gData$markers),GWA.result,allele.frequency=all.allele.frequencies)

      if (ncol(gData$genes)!=0) {
        old.names           <- colnames(GWA.result)
        GWA.result          <- cbind(GWA.result,gData$map$gene1,gData$map$gene2)
        colnames(GWA.result)<- c(old.names,"gene1","gene2")
      }

      ########################


      if (max.score==1) {
        GWA.result$effect   <- 0.5 * GWA.result$effect
      }

      ##################

      write.table(GWA.result,file=output.file,quote=F,row.names=F,col.names=T,sep=",")
      rm(X.temp)
    }


    N.effective <- sum(!is.na(GWA.result$pvalue))
    analyzed.snps <- which(!is.na(GWA.result$pvalue))

    ##########################################
    #CALCULATE A SIGNIFICANCE THRESHOLD

    # When BT is 1,2,3 or 4, determine the LOD-threshold
    if (BT==1) {LOD.thr <- -log10(alpha/N.effective) }# -log10(alpha/gData$N)}
    #if (BT==2) {LOD.thr <- -log10(1/gData$N)}
    if (BT==3) {LOD.thr <- sort(-log10(na.omit(GWA.result$pvalue)),decreasing=T)[K]}

    if (BT==4) {
      cut.off <- 0.995
      number.of.nonmissing  <- aggregate(gData$pheno[,tr.n],by=list(ordered(gData$pheno$genotype)),FUN=function(x){sum(!is.na(x))})[match(gData$plant.names,sort(gData$plant.names)),2]
      number.of.nonmissing[number.of.nonmissing>0] <- 1
      ind.indices <- rep((1:gData$n)[number.of.nonmissing>0],times=number.of.nonmissing[number.of.nonmissing>0])
      SIGMA       <- varcomp.values[1,] * gData$kinship[ind.indices,ind.indices] + varcomp.values[2,] * diag(sum(number.of.nonmissing))
      COR         <- cov2cor(SIGMA)
      INV.COR     <- GINV(COR)
      Keff      <- 0
      n.block   <- 0
      b.size    <- 10*sum(number.of.nonmissing)
      for (CHR in 1:gData$nchr) {
        blocks  <- DefineBlocks(which(gData$map$chromosome==CHR),block.size=b.size)
        n.block <- n.block  +  length(blocks)
        for (b in 1:length(blocks)) {
          marker.frame <- gData$markers[blocks[[b]],]
          Keff <- Keff + GaoCorrection(marker.frame=marker.frame,number.of.replicates=number.of.nonmissing,inv.cor.matrix=INV.COR,cut.off=cut.off)
        }
      }
      LOD.thr <- -log10(alpha/Keff)
    }

    #######################################################################
    #GIVEN THE SIGNIFICANCE THRESHOLD FROM THE PREVIOUS STEP, SELECT THE SNPs WHOSE LOD-SCORE IS ABOVE THIS THRESHOLD

    snp.selection         <- (-log10(ReplaceNaByOne(GWA.result$pvalue)) >= LOD.thr)
    snp.selection[is.na(snp.selection)] <- F

    ########################################################################

    if (sum(snp.selection,na.rm=T) > 0) {
      #
      selection.numbers <- which(snp.selection)
      snp.selection.extra   <- (1:gData$N %in% selection.numbers)

      if (sizeIncludedRegion > 0) {
        if (extra.output.file & nrow(gData$genes)>0) {
          for (mg in which(snp.selection)) {
            #selection.numbers <- c(selection.numbers,GetSNPsInRegion(gData,mg,sizeIncludedRegion))
            selection.numbers <- c(selection.numbers,
              GetSNPsInRegionWithSufficientLD(gData,snp.number=mg,
                region.size=sizeIncludedRegion,
                minR2=minR2))
          }
          selection.numbers           <- sort(unique(selection.numbers))
          snp.status                  <- rep(paste("within",as.character(sizeIncludedRegion/1000),"kb of a significant snp"),length(selection.numbers))
          significant.among.selection <- which(selection.numbers %in% which(snp.selection))
          snp.status[significant.among.selection] <- "significant snp"

          snp.selection            <- (1:gData$N %in% selection.numbers)

          selection.numbers.extra  <- which(snp.selection.extra)
          snp.selection.extra      <- (1:gData$N %in% selection.numbers.extra)

          snp.sig.extra               <- data.frame(marker=row.names(gData$markers)[snp.selection.extra],
            chromosome=gData$map$chromosome[snp.selection.extra],
            position=gData$map$position[snp.selection.extra],
            pvalue=GWA.result$pvalue[snp.selection.extra])

          snp.sig.extra               <- data.frame(snp.sig.extra,gene1=as.character(gData$map$gene1[snp.selection.extra]),
            gene2=as.character(gData$map$gene2[snp.selection.extra]),
            other.genes=as.character(rep("",length(selection.numbers.extra))))

          snp.sig.extra$other.genes   <- as.character(snp.sig.extra$other.genes )

          for (mg in which(snp.selection.extra)) {
            mg.index <- which(which(snp.selection.extra)==mg)
            extra.snps <- GetSNPsInRegionWithSufficientLD(gData,snp.number=mg,region.size=sizeIncludedRegion,minR2=minR2)
            if (length(extra.snps) > 0) {
              new.genes <- na.omit(as.character(sort(unique(gData$map$gene1[extra.snps],unique(gData$map$gene2[extra.snps])))))
              snp.sig.extra$other.genes[mg.index] <- paste(setdiff(new.genes,na.omit(c(as.character(snp.sig.extra$gene1[mg.index]),as.character(snp.sig.extra$gene2[mg.index])))),collapse=',')
            }
          }
          snp.output.file.extra       <- paste("results/","significant.snps.with.extra.genes_",trait,suffix,".csv",sep="")
          write.table(snp.sig.extra,file=snp.output.file.extra,quote=FALSE,row.names=FALSE,sep=",")

        } else {
          for (mg in which(snp.selection)) {
            selection.numbers <- c(selection.numbers,GetSNPsInRegionWithSufficientLD(gData,snp.number=mg,region.size=sizeIncludedRegion,minR2=minR2))
          }
          selection.numbers           <- sort(unique(selection.numbers))
          snp.status                  <- rep(paste("within",as.character(sizeIncludedRegion/1000),"kb of a significant snp"),length(selection.numbers))
          significant.among.selection <- which(selection.numbers %in% which(snp.selection))
          snp.status[significant.among.selection] <- "significant snp"
          snp.selection   <- (1:gData$N %in% selection.numbers)
        }
      } else {
        snp.status                  <- rep("significant snp",length(selection.numbers))
      }


      snp.sig               <- data.frame(marker=row.names(gData$markers)[snp.selection],
        chromosome=gData$map$chromosome[snp.selection],
        position=gData$map$position[snp.selection],
        pvalue=GWA.result$pvalue[snp.selection],
        snp.status=snp.status)


      if (nrow(gData$genes)>0) {

        snp.sig               <- data.frame(snp.sig,gene1=gData$map$gene1[snp.selection],gene2=gData$map$gene2[snp.selection])

        if (gData$nchr==5) { # if the species is arabisopsis ...
          gene1http <- as.character(snp.sig$gene1)
          gene2http <- as.character(snp.sig$gene2)
          gene1http[!is.na(gene1http)] <- paste("=HYPERLINK(\"http://arabidopsis.org/servlets/TairObject?type=locus&name=",gene1http[!is.na(gene1http)],"\")",sep="")
          gene2http[!is.na(gene2http)] <- paste("=HYPERLINK(\"http://arabidopsis.org/servlets/TairObject?type=locus&name=",gene2http[!is.na(gene2http)],"\")",sep="")
          snp.sig               <- data.frame(snp.sig,gene1_http=gene1http,gene2_http=gene2http)
          snp.sig$gene1_http<- as.character(snp.sig$gene1_http)
          snp.sig$gene2_http<- as.character(snp.sig$gene2_http)
        }

      }

      snp.sig$marker        <- as.character(snp.sig$marker)
      no.significant.snps   <- FALSE

      allele.frequencies    <- all.allele.frequencies[snp.selection]
      #apply(gData$markers[snp.sig$marker,unique(gData$pheno$genotype[!is.na(gData$pheno[,tr.n])])],1,mean)

      #if (max.score==2) {allele.frequencies <- allele.frequencies / 2}


      if (GLSMethod==1) {
        effect.sizes                  <- GWA.result$effect[snp.selection]

        variance.explained.by.snp     <- 4 * effect.sizes^2 * all.allele.variances[snp.selection]
        proportion.var.snp.fixed <- variance.explained.by.snp / var(as.numeric(na.omit(gData$pheno[,tr.n])))

        snp.sig             <- cbind(snp.sig,data.frame(reference.allele.freq=allele.frequencies,
          effect.size=effect.sizes,
          proportion.var.snp.fixed=proportion.var.snp.fixed))
        #perc.of.phenetic.var=explained.phenotypic.variance))
      }

      if (GLSMethod==2) {
        effect.sizes                  <- GWA.result$effect[snp.selection]
        #variance.explained.by.snp     <- 4 * effect.sizes^2 * allele.frequencies * (1 - allele.frequencies)
        #explained.phenotypic.variance <- 100 * variance.explained.by.snp / var(as.numeric(na.omit(gData$pheno[,tr.n])))
        #
        snp.sig             <- cbind(snp.sig,data.frame(reference.allele.freq=allele.frequencies,
          effect.size=effect.sizes,
          proportion.var.LRT=GWA.result$R_LR_2[snp.selection]))
        #                             perc.of.phenetic.var=explained.phenotypic.variance))
      }

      if (GLSMethod %in% 3:4) {

        effect.sizes          <- GWA.result$effect[snp.selection]
        #variance.explained.by.snp       <- 4 * effect.sizes^2 * allele.frequencies * (1 - allele.frequencies)
        #explained.phenotypic.variance   <- 100* variance.explained.by.snp / var(as.numeric(na.omit(gData$pheno[,tr.n])))

        variance.explained.by.snp     <- 4 * effect.sizes^2 * all.allele.variances[snp.selection]
        proportion.var.snp.fixed      <- variance.explained.by.snp / var(as.numeric(na.omit(gData$pheno[,tr.n])))

        snp.sig             <- cbind(snp.sig,data.frame(reference.allele.freq=allele.frequencies,
          effect.size=effect.sizes,
          proportion.var.LRT=GWA.result$R_LR_2[snp.selection],
          proportion.var.snp.fixed=proportion.var.snp.fixed))
        #                             perc.of.phenetic.var=explained.phenotypic.variance))
      }


      if (var.explained.snp.random) {

        require(asreml)
        K.main              <- GRM(t(as.matrix(gData$markers)))
        Kinv.asreml         <- MakeKinshipAsreml(K.main)

        # the following lines are in fact redundant if asreml was already used before to estimate the variance components.
        if (sum(cov.cols)!=0) {  # reml.formula is the formula for the fixed part of the model
          reml.formula        <- as.formula(paste(paste(names(gData$pheno)[tr.n],"~"),paste(names(gData$pheno)[cov.cols],collapse="+")))
        } else {
          reml.formula        <- as.formula(paste(names(gData$pheno)[tr.n],"~ 1"))
        }


        random.result <- rep(0,nrow(snp.sig))

        for (i in 1:nrow(snp.sig)) {
          unstandardized.snp  <- gData$markers[snp.sig$marker[i],]
          standardized.snp    <- as.numeric(scale(as.numeric(unstandardized.snp)))
          names(standardized.snp) <- names(unstandardized.snp)

          reml.obj.i <- asreml(fixed= reml.formula,random = as.formula(" ~ giv(genotype,var=T) + snp"),
            data=data.frame(gData$pheno[!is.na(gData$pheno[,tr.n]),],
              snp=standardized.snp[gData$pheno$genotype[!is.na(gData$pheno[,tr.n])]]),
            na.method.X="omit",ginverse =list(genotype=Kinv.asreml))

          smv <- summary(reml.obj.i)$varcomp$component
          random.result[i] <- (smv[2] / sum(smv))
        }
        snp.sig             <- data.frame(snp.sig,proportion.var.snp.random=random.result)
        colnames(snp.sig)[ncol(snp.sig)] <- 'proportion.var.snp.random'
      }

      #############

      # restore for arabidopsis ?
      #accessions.with.the.rare.allele <- rep("",nrow(snp.sig))
      #for (mg in 1:nrow(snp.sig)) {
      #  if (snp.sig$reference.allele.freq[mg] < 0.5) {rare.allele<-1} else {rare.allele<-0}
      #  temp.obj <- gData$markers[snp.sig$marker[mg],unique(gData$pheno$genotype[!is.na(gData$pheno[,tr.n])])]
      #  temp.names <- names(temp.obj)
      #  accessions.with.the.rare.allele[mg] <- paste(temp.names[temp.obj==rare.allele],collapse=",")
      #}
      #snp.sig <- cbind(snp.sig,data.frame(accessions.with.the.rare.allele=accessions.with.the.rare.allele))

      # creating this file
      snp.output.file       <- paste("results/","significant.snps.",trait,suffix,".csv",sep="")

      write.table(snp.sig,file=snp.output.file,quote=FALSE,row.names=FALSE,sep=",")

      if (make.all.trait.significant.snp.file) {
        snp.sig <- data.frame(snp.sig[,1:(ncol(snp.sig)-1)],trait.name=rep(trait,nrow(snp.sig)),snp.sig[,ncol(snp.sig)])
        write.table(snp.sig,file=paste("results/",all.trait.significant.snp.file,sep=""),append=T,quote=FALSE,row.names=FALSE,col.names=F,sep=",")
      }


    } else {

      no.significant.snps         <- TRUE

    }

    if (make.all.trait.h2.file) {
      all.traits$h2[which(trait.numbers==tr.n)] <- (varcomp.values[1,1])/(sum(varcomp.values[,1])) # which(names(gData$pheno)==trait)-min(trait.numbers)+1
      all.traits$inflation[which(trait.numbers==tr.n)] <- inflation.factor
      if (kinship.type %in% 3:4) {all.traits$n.markers.in.K[which(trait.numbers==tr.n)] <- testK$m.optimal}
    }

    if (make.one.file.with.all.lod.scores) {
      all.lod.scores.frame[,which(names(all.lod.scores.frame)==trait)] <- -log10(ReplaceNaByOne(GWA.result$pvalue))
    }

    if (make.one.file.with.all.snp.effects) {
      if (GLSMethod==1) {all.snp.effects.frame[,which(names(all.snp.effects.frame)==trait)] <- GWA.result$effect}
      if (GLSMethod %in% 2:4) {all.snp.effects.frame[,which(names(all.snp.effects.frame)==trait)] <- GWA.result$effect}
    }

    ############################################################################

    # Create the summary file:

    summ.file    <- paste("results/","summary.",trait,suffix,".txt",sep="")


    cat("R-image used: ",r.image.name,"\n","\n",file=summ.file)     # file is created and first line is written; APPEND=FALSE (default)
    cat("Trait: ",trait,"\n","\n",file=summ.file,append=T)
    cat("Analysis started on: ",start.date,"\n",file=summ.file,append=TRUE)
    cat("Analysis finished on: ",date(),"\n","\n",file=summ.file,append=TRUE)

    cat("Data are available for",gData$N,"SNPs", "\n",file=summ.file,append=TRUE)
    if (MAF > 0) {cat(gData$N-N.effective,"of them were not analyzed because their minor allele frequency is below",MAF,"\n","\n",file=summ.file,append=TRUE)}


    if (GLSMethod %in% c(1,3)) {
      cat("Mixed model with only polygenic effects, and no marker effects:","\n",file=summ.file,append=TRUE)
      cat("Genetic variance: ",varcomp.values[1,1],"\t","standard error: ",varcomp.std[1],"\n",file=summ.file,append=TRUE)
      cat("Residual variance: ",varcomp.values[nrow(varcomp.values),1],"\t","standard error: ",varcomp.std[nrow(varcomp.values)],"\n\n",file=summ.file,append=TRUE)
      #
      cat("File containing the p-values of all snps:  ",output.file,"\n","\n",append=TRUE,file=summ.file,sep="")
    }
    #

    if (BT %in% 1:4) {
      cat("LOD-threshold: ",LOD.thr,"\n",append=TRUE,file=summ.file)
      if (BT==4) {cat("Number of effective tests: ",Keff,"\n",append=TRUE,file=summ.file)}
      #
      if (!no.significant.snps) {
        cat("File containing the p-values of the selected snps: ",snp.output.file,"\n",append=TRUE,file=summ.file)
        cat("Number of selected snps =",nrow(snp.sig),"\n",append=TRUE,file=summ.file)
        cat("Smallest p-value among the selected snps:",min(snp.sig$pvalue),"\n",append=TRUE,file=summ.file)
        cat("Largest  p-value among the selected snps:",max(snp.sig$pvalue),"(LOD-score:",-log10(max(snp.sig$pvalue)),")","\n",append=TRUE,file=summ.file)
      } else {
        cat("No significant snps found.","\n",append=TRUE,file=summ.file)
      }
      if (genomic.control) {cat("\n","Genomic control correction was applied","\n",append=TRUE,file=summ.file)}
      if (!genomic.control) {cat("\n","No Genomic control correction","\n",append=TRUE,file=summ.file)}
      cat("Genomic control inflation-factor = ",GC[[2]],"\n","\n",append=TRUE,file=summ.file)
    }


    if ((sum(cov.cols)!=0)  &  reml.algo == "asreml") {
      cat(paste("pvalues for the covariables",paste(names(gData$pheno)[cov.cols],collapse=","),"in the mixed model without snps:"),"\n",append=TRUE,file=summ.file)
      cat((wald(reml.obj))[1+(1:length(cov.cols)),4],"\n","\n",append=TRUE,file=summ.file)
    }

    ##################################################################################################

    # Make a plot of the LOD-profile: (pdf and jpeg)
    LOD.thr.plot  <- 0
    if (BT %in% 1:4) {LOD.thr.plot <- LOD.thr}
    if (BT == 5) {if (!no.significant.snps) {LOD.thr.plot <- min(-log10(GWA.result$pvalue[snp.selection]))}}
    #
    x.effects=integer(0)
    effects.size=numeric(0)

    if ("real.effects" %in% names(gData)) {
      if (nrow(gData$real.effects$locations)>0) {
        x.effects=gData$real.effects$locations[,tr.n-1]
        effects.size=gData$real.effects$sizes[,tr.n-1]
      }
    }

    if (no.significant.snps) {
      x.sig <- integer(0)
    } else {
      x.sig <- which(snp.selection)[which(snp.sig$snp.status=="significant snp")]
    }


    MakeLodPlotWithChromosomeColors(xvalues=gData$map$cum.position,yvalues=-log10(ReplaceNaByOne(GWA.result$pvalue)),plot.title=trait,gData=gData,
      col.palette = rep(c("royalblue","maroon"),50)[1:gData$nchr],file.name=paste("plots/",trait,suffix,".jpeg",sep=""),
      x.sig=x.sig,chr.boundaries=cumsum(gData$chr.lengths.bp)[-1],y.thr=LOD.thr.plot,
      x.effects=x.effects,effects.size=effects.size,plot.type="d")

    if (!jpeg.only) {
      MakeLodPlotWithChromosomeColors(xvalues=gData$map$cum.position,yvalues=-log10(GWA.result$pvalue),plot.title=trait,gData=gData,
        col.palette = rep(c("royalblue","maroon"),50)[1:gData$nchr],
        file.name=paste("plots/",trait,suffix,".pdf",sep=""),x.sig=x.sig,
        chr.boundaries=cumsum(gData$chr.lengths.bp)[-1],y.thr=LOD.thr.plot,x.effects=x.effects,
        effects.size=effects.size,jpeg.plot=FALSE,plot.type="d")
    }

    QQplotPvalues2(pvalues=GWA.result$pvalue,main.title=trait,file.name=paste0("plots/QQplot_",trait,suffix,".png"))


  }   # end for (tr.n in trait.numbers)


  ########################################################################################################################
  ########################################################################################################################

  if (make.all.trait.h2.file) {
    write.csv(all.traits,file=paste("results/",all.trait.h2.file,sep=""),row.names=F)
  }

  if (make.one.file.with.all.lod.scores) {
    write.csv(all.lod.scores.frame,file=paste("results/",file.with.all.lod.scores,sep=""),quote=F)
  }

  if (make.one.file.with.all.snp.effects) {
    write.csv(all.snp.effects.frame,file=paste("results/",file.with.all.snp.effects,sep=""),quote=F)
  }

  #if (covariables) {
  #  # After the loop over traits
  #  gData$pheno <- temp.pheno
  #}
}
