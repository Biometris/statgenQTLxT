add_Phenotypic_Data_To_Gwas_Object <- function(csv.file.name,r.image.name,new.r.image.name,data.path=setwd(),which.columns.as.factor=integer(0)) {
# OBJECTIVE :
#' @param
#' @param
# INPUT :
#' @param csv.file.name : name of the csv file containing the data (first column-name: genotype). Contains no row-names
#' @param r.image.name,new.r.image.name :
#' @param data.path : the path where the function will look for the above files
#' @param
#' @param which.columns.as.factor : column-numbers of the variables that are to be imported as factor
# OUTPUT :
# *
# *
# *
#NB. side effect : working dit. is changed

#gwas.obj=GWAS.obj;csv.file.name=csv.file.name;add.var.means=F;make.pheno.image=F;pheno.image.name=""
  #setwd(data.path)
  load(r.image.name)
  gwas.obj <- AddPhenoData(gwas.obj=gwas.obj,csv.file.name=csv.file.name,add.var.means=F,make.pheno.image=F,pheno.image.name="",
                           which.columns.as.factor=which.columns.as.factor)
  for (k in which.columns.as.factor) {
      if ("NA" %in% levels(gwas.obj$pheno[,k])) {
      gwas.obj$pheno[which(is.na(gwas.obj$pheno[,k])),k] <- "NA"
      is.na(gwas.obj$pheno[,k])[gwas.obj$pheno[,k]=="NA"] <- TRUE
      suppressWarnings(gwas.obj$pheno[,k] <- factor(gwas.obj$pheno[,k],exclude="NA"))
    }
  }
  save(gwas.obj,file=new.r.image.name)
}

GetSNPsInRegion <- function(gwas.obj,snp.number,region.size=5000){
# OBJECTIVE : GIVE THE SNP-NUMBERS of the snps that are within region.size base-pairs (on either side) of the target snp with number 'snp.number'
# both snp.number and the output refer to row-numbers in gwas.obj$markers; not to base pair positions
  crit1 <- (abs(gwas.obj$map$position[snp.number] - gwas.obj$map$position) <= region.size)
  crit2 <- (gwas.obj$map$chromosome==gwas.obj$map$chromosome[snp.number])
  return(sort(setdiff(which(crit1 & crit2),snp.number)))
}

GetSNPsInRegionWithSufficientLD <- function(gwas.obj,snp.number,region.size=5000,min.r2=0.5,snps.as.columns=F){
# OBJECTIVE : GIVE THE SNP-NUMBERS of the snps that are within region.size base-pairs (on either side) of the target snp with number 'snp.number'
# both snp.number and the output refer to row-numbers in gwas.obj$markers; not to base pair positions
#
# Refinement (compared to GetSNPsInRegion) : only those snps that are in sufficient LD with the target snp (snp.number) are returned
# (where LD is measured in terms of r^2)
#
# snps.as.columns : put to TRUE when in GWAS.obj$markers the markers are in the columns ('n x p format')
#
# gwas.obj=GWAS.obj;snp.number=10000;region.size=50000;min.r2=0.3 ; snps.as.columns=T
  crit1 <- (abs(gwas.obj$map$position[snp.number] - gwas.obj$map$position) <= region.size)
  crit2 <- (gwas.obj$map$chromosome==gwas.obj$map$chromosome[snp.number])
  candidate.snps <- sort(setdiff(which(crit1 & crit2),snp.number))
  if (snps.as.columns) {
    r2.matrix <- cor((as.matrix(gwas.obj$markers[,candidate.snps])),(as.matrix(gwas.obj$markers[,snp.number])))^2
    candidate.snps.names <- names(which(r2.matrix[,1]> min.r2))
    return(sort(  which(colnames(gwas.obj$markers) %in% candidate.snps.names)))
  } else {
    r2.matrix <- cor(t(as.matrix(gwas.obj$markers[candidate.snps,])),t(as.matrix(gwas.obj$markers[snp.number,])))^2
    candidate.snps.names <- names(which(r2.matrix[,1]> min.r2))
    return(sort(  which(row.names(gwas.obj$markers) %in% candidate.snps.names)))
  }
}
# GetSNPsInRegionWithSufficientLD(GWAS.obj,snp.number=10000,region.size=100000,min.r2=0.5,snps.as.columns=T)

MakeSubsampleFrame <- function(phenotype.name,subsample.vector,population.vector,genotype.vector,phenotype.vector) {
  population.vector <- as.integer(population.vector)
  subsample.vector  <- as.integer(subsample.vector)
  not.sampled <- setdiff(population.vector,subsample.vector)
  #
  data.sub                 <- data.frame(genotype.vector,phenotype.vector)
  names(data.sub)          <- c("genotype",phenotype.name)
  data.sub[not.sampled,2]  <- NA
return(data.sub)
}

PrepareScanGls <- function(gwas.obj,tr.n,cov.cols=integer(0),data.path=getwd(),suffix="") {
#browser()
  if (sum(cov.cols)!=0) {
    reml.formula        <- as.formula(paste(paste(names(gwas.obj$pheno)[tr.n],"~"),
                                      paste(names(gwas.obj$pheno)[cov.cols],collapse="+")))
  } else {
    reml.formula        <- as.formula(paste(names(gwas.obj$pheno)[tr.n],"~ 1"))
  }
  #reml.obj              <- asreml(maxiter = 25,fixed= reml.formula,data=gwas.obj$pheno, random = ~ giv(genotype,var=T),
  #                                na.method.X="omit",ginverse = list(genotype=gwas.obj$kinship.asreml))
  reml.obj              <- asreml(maxiter = 25,fixed= reml.formula,data=gwas.obj$pheno, random = ~ giv(genotype,var=T),na.method.X="omit",ginverse = list(genotype=gwas.obj$kinship.asreml))
  #asreml(fixed= reml.formula,random = ~ giv(genotype,var=T),data=gwas.obj$pheno,sparse= ~1,asreml.control(maxiter = 25,ginverse = list(genotype = gwas.obj$kinship.asreml)),na.method.X="omit")
  varcomp.values      <- data.frame(var.comp.values=summary(reml.obj)$varcomp$component)
  varcomp.file        <- paste(data.path,"output/",trait,".","varcomp",".csv",sep="")
  MakeVarcompFile(var.comp.values=varcomp.values,file.name=varcomp.file)
return(list(varcomp.values=varcomp.values,varcomp.file=varcomp.file,reml.obj=reml.obj,reml.formula=reml.formula))
}
# temp.obj        <- PrepareScanGls(gwas.obj=GWAS.obj,tr.n=tr.n,cov.cols=cov.cols,data.path=data.path,suffix=suffix)

MakePlinkFiles <- function(gwas.obj,file.name="gwas",make.bed=T) {    # data.path=getwd(),plink.path=data.path # plink.path : path of the plink executable/binary # data.path  : path where the files are to be created
# file.name  : file-name (without tped,bim etc extension)
# N.B. the plink executable/binary should be in the current working directory
# OUTPUT : ...
# N.B. standard ped/map files can be obtained from tped/tfam using plink : plink --tfile oldname --recode newname       # or  :--out newname ??
# , where old/new name are without .tped etc extension
#
# from tfam/tped to bim/bed/fam :
# plink --tfile LFN349acc_001 --make-bed --out LFN349acc_001

  # TO DO : WINDOWS/LINUX
  #stopifnot(file.exists("plink.exe"))

  block.list <- DefineBlocks(1:gwas.obj$N,5000)
  tped.name  <- paste(file.name,".tped",sep="") # data.path
  cat("",file=tped.name)
  if (is.integer(gwas.obj$map$position)) {position.in.cm <- F} else {position.in.cm <- T}
  for (b in 1:length(block.list)) {
    block <- block.list[[b]]
    if (max(gwas.obj$markers)==1) {
      #snp.out  <- gwas.obj$markers[block,rep(match(gwas.obj$pheno$genotype,gwas.obj$plant.names),each=2)] + 1
      # 15 may 2015:
      snp.out  <- gwas.obj$markers[block,rep(gwas.obj$plant.names,each=2)] + 1
    } else {
      snp.out  <- gwas.obj$markers[block,rep(gwas.obj$plant.names,each=2)]
      for (marker in block) {
        het.ind  <- which(gwas.obj$markers[marker,]==1)
        snp.out[which(block==marker),het.ind*2] <- 2
        snp.out[which(block==marker),(het.ind*2 - 1) ] <- 0
      }
    }
    if (position.in.cm) {
      output.frame    <- data.frame(gwas.obj$map$chromosome[block],row.names(gwas.obj$markers)[block],gwas.obj$map$position[block],matrix(0,length(block),1),snp.out)
    } else {
      output.frame    <- data.frame(gwas.obj$map$chromosome[block],row.names(gwas.obj$markers)[block],matrix(0,length(block),1),gwas.obj$map$position[block],snp.out)
    }
    write.table(output.frame, file = tped.name, append = T, quote = FALSE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)
  }

  if (make.bed) {
    gwas.obj$pheno[,2] <- rnorm(nrow(gwas.obj$pheno),mean=20)
    MakeTfamFile(gwas.obj,accessions.as.families=T,file.name=file.name)
    #system(paste("plink --tfile",file.name,"--make-bed"))
    command <- paste("plink --tfile",file.name,"--make-bed","--out",file.name)
    if (position.in.cm) {command <- paste(command,"--cm")}
    system(command)
  }
}

# MakeTfamFile(GWAS.obj,trait.nr  = cov.cols,file.name=covar.file,file.type=1,delimiter="\t")
# MakeTfamFile <- function(gwas.obj,trait.nrs  = 2,accessions.as.families=F,file.names= names(gwas.obj$pheno)[trait.nrs]) {
MakeTfamFile <- function(gwas.obj,trait.nr  = 2,accessions.as.families=T,file.name=gwas.obj$description,file.type=0,delimiter=" ") {

# N.B. only works for arabidopsis data; NOT for WTCCC data
# TO DO : adapt to out-breeders

# file.name : file-name, without tped extension
# file.type : 0 = 6 column tped format;
#             1 = 3 column txt ("alternate phenotype") format
  #pn  <- length(plant_names)
  n.ind <- nrow(gwas.obj$pheno)
  if (accessions.as.families) {
    fam.vector <- gwas.obj$pheno$genotype
  } else {
    fam.vector <- 1:n.ind
  }
  trait.values <- matrix(-9,n.ind,length(trait.nr))
  for (tn in trait.nr) {
    trait.values[!is.na(gwas.obj$pheno[,tn]),which(tn==trait.nr)] <- gwas.obj$pheno[!is.na(gwas.obj$pheno[,tn]),tn]
  }
  trait.values <- as.data.frame(trait.values)
  names(trait.values) <- names(gwas.obj$pheno[trait.nr])
  if (file.type==0) {
    #write.table(data.frame(fam.vector,1:n.ind,rep(0,n.ind),rep(0,n.ind),rep(1,n.ind),trait.values), file = paste(file.name,".tfam",sep=""), quote = FALSE, sep = " ",row.names = FALSE,col.names = FALSE)
    write.table(data.frame(fam.vector,row.names(gwas.obj$pheno),rep(0,n.ind),rep(0,n.ind),rep(1,n.ind),trait.values),
                file = paste(file.name,".tfam",sep=""), quote = FALSE, sep = delimiter,row.names = FALSE,col.names = FALSE)
  } else {
    #write.table(data.frame(fam.vector,1:n.ind,trait.values), file = paste(file.name,".txt",sep=""), quote = FALSE, sep = " ",row.names = FALSE,col.names = FALSE)
    write.table(data.frame(fam.vector,row.names(gwas.obj$pheno),trait.values), file = paste(file.name,".txt",sep=""), quote = FALSE, sep = delimiter,row.names = FALSE,col.names = FALSE)
  }
}
#MakeTfamFile(GWAS.obj)
# to make a Fast-LM covariate - file :
#MakeTfamFile(GWAS.obj,trait.nr  = 3:5,file.name=trait,file.type=1,delimiter="\t")

AddPhenoData2  <- function(gwas.obj,csv.file.name,add.var.means=FALSE,mean.cols=0,make.pheno.image=F,pheno.image.name="pheno.RData",add.normal.transform=F) {
# assumes that the working directory is correct, and that the functions script is loaded
# input  : the "gwas-object" gwas.obj and the phenotypic csv-file csv.file.name
# output : the same gwas.obj, the element gwas.obj$pheno containing the data from the csv file
#           if make.pheno.image=TRUE, also an R-image containing the phenotypic data is created
#           (for consistency with earlier version this data-frame is then called data1)

  if (!add.var.means) {mean.cols<-0}
  if (sum(mean.cols)==0) {add.var.means<-FALSE}
  #
  data0                   <- read.table(file=csv.file.name,sep=",",na.strings="",header=T)
  data0                   <- read.table(file=csv.file.name,sep=",",na.strings="",header=T,colClasses=c("character",rep("numeric",ncol(data0)-1))) # ,colClasses=rep("numeric",12)
  # complete the genotype-names:
  for (i in 2:nrow(data0)) {if (is.na(data0[i,1])) {data0[i,1] <- data0[i-1,1]}}
  names(data0)[1]         <- "genotype"
  # convert the accession names from factor to character
  data0$genotype          <- as.character(data0$genotype)
  data1                   <- data.frame(data0, stringsAsFactors = F)
  # remove the accessions which are not in plant.names
  data1                   <- data1[setdiff(1:nrow(data1),which(data1$genotype %in% setdiff(unique(data1$genotype),gwas.obj$plant.names))),]
  # check the number of replicates in the remaining genotypes
  #rep.vec                 <- tabulate(factor(data1$genotype))
  rep.vec                 <- as.numeric(table(data1$genotype)[match(unique(data1$genotype),sort(unique(data1$genotype)))])


  # Add letters a,b,c... to the row-names
  ##!##
  row.names(data1)[1:(rep.vec[1])]<- paste(data1$genotype[1:(rep.vec[1])], letters[1:(rep.vec[1])],sep="")
  for (i in 2:length(rep.vec)) {row.names(data1)[sum(rep.vec[1:(i-1)])+ 1:(rep.vec[i])]<- paste(data1$genotype[sum(rep.vec[1:(i-1)])+ 1:(rep.vec[i])], letters[1:(rep.vec[i])],sep="")}
  # Now the row names are the accession names, where each accession name is replicated, with lower case letters as extension,
  # e.g. CS76113a,CS76113b,CS76113c

  data2                   <- data1
  # Now add extra accessions, i.e. the ones that do not occur in data1/data2, but DO occur in plant.names(2?)
  n.extra <- length(setdiff(gwas.obj$plant.names,unique(data1$genotype)))
  if (n.extra>0) {
    extra.acc               <- data.frame(matrix(NA,n.extra,ncol(data2)))
    extra.acc[,1]           <- setdiff(gwas.obj$plant.names,unique(data2$genotype))
    row.names(extra.acc)    <- paste(extra.acc[,1],"a",sep="")
    names(extra.acc)        <- names(data2)
    data2                   <- rbind(data2,extra.acc)
    rep.vec                 <- c(rep.vec,rep(1,n.extra))
  }


  rep.vec2     <- rep.vec[match(gwas.obj$plant.names,unique(data2$genotype))]
  plant.names2 <- rep(gwas.obj$plant.names,times=rep.vec2)
  plant.names2[1:(rep.vec2[1])] <- paste(plant.names2[1:(rep.vec2[1])], letters[1:(rep.vec2[1])],sep="")

  for (i in 2:length(rep.vec2)) {
    plant.names2[sum(rep.vec2[1:(i-1)])+ 1:(rep.vec2[i])] <-
    paste(plant.names2[sum(rep.vec2[1:(i-1)])+ 1:(rep.vec2[i])], letters[1:(rep.vec2[i])],sep="")
  }

  data3   <- data2[match(plant.names2,row.names(data2)),]

  data1   <- data3
  rep.vec <- rep.vec2


  if (add.var.means) {
    NC    <- ncol(data1)
    data1 <- AddMeans(input.frame=data1,col.select=mean.cols)   # the averages of the columns given in col.select are now added as extra columns
    if (add.normal.transform) {
      for (i in c(mean.cols,NC+1:length(mean.cols))) {
        data1 <- cbind(data1,qqnorm(data1[,i],plot.it=F)$x)
        names(data1)[ncol(data1)] <- paste(names(data1)[i],"_transformed",sep="")
      }
    }
  }


  if (make.pheno.image) {save(data1,rep.vec,file=pheno.image.name)}


  gwas.obj$pheno  <- data1
  # old code :
  #save(gwas.obj,file=new.gwas.image)
  #assign(gwas.obj$pheno, value=data1, envir = .GlobalEnv)

return(gwas.obj)
}




# function found on the internet: (author?)
is.installed <- function(mypkg) {return(is.element(mypkg, installed.packages()[,1]))}

# Taken from : http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
MatrixRoot <- function(x) { # assumes that x is symmetric
  x.eig <- eigen(x,symmetric=TRUE)
  x.sqrt <- x.eig$vectors %*% diag(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  return(x.sqrt)
}

# Taken from : http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
denman.beavers <- function(mat,maxit=50) {
  stopifnot(nrow(mat) == ncol(mat))
  niter <- 0
  y <- mat
  z <- diag(rep(1,nrow(mat)))
  for (niter in 1:maxit) {
    y.temp <- 0.5*(y+solve(z))
    z <- 0.5*(z+solve(y))
    y <- y.temp
  }
  return(list(sqrt=y,sqrt.inv=z))
}

MatrixRoot2 <- function(x) {
  ev <- eigen(x)
  V <- ev$vectors
  return(V %*% diag(sqrt(ev$values)) %*% t(V))
}


asv     <- function(x) {return(summary(as.vector(x)))}


GINV    <- function(M) {
  svdM    <- svd(M)
  return(svdM$v %*%diag(1/svdM$d)%*% t(svdM$u))
}

GRM     <- function(X) {
# X = n x p matrix of marker scores
  X <- as.matrix(X)
  N <- ncol(X)
  require(pls)
  Y <- stdize(X)
  K <- Y %*% t(Y) / N
  K <- K / KinshipTransform(K)
return(K)
}


IBS     <- function(X,normalization=nrow(X)) {
                X       <- as.matrix(X)
                Y       <- X
                Y[X==0] <- 1
                Y[X==1] <- 0
return((t(X) %*% X + t(Y) %*% Y)/normalization)
}

IBS3     <- function(X,normalization=4*nrow(X)) {
                X       <- as.matrix(X)
                X       <- replace(X, list=1:2, values=2:1)
                Y       <- X
                Y[X%in% 1:2] <- 0
                Y[X==0] <- 2
                #Y       <- replace(X, list=0:2, values=c(2,0,0))
                #Z       <- replace(X, list=0:2, values=c(0,2,0))
                #X       <- replace(X, list=0:2, values=c(0,0,2))
                Z       <- X
                Z[X%in% c(0,2)] <- 0
                Z[X==1] <- 2
                X[X==1] <- 0
                return((t(X) %*% X + t(Y) %*% Y + t(Z) %*% Z)/normalization)
}
# IBS(GWAS_obj$markers[1:10,])[1:16,1:10]- IBS3(GWAS_obj$markers[1:10,])[1:16,1:10]
# FIND THE MISTAKE !!

IBS2    <- function(X1,X2,normalization= max(nrow(X1),nrow(X2),ncol(X1),ncol(X2))) {
                X1       <- as.matrix(X1)
                Y1       <- X1
                Y1[X1==0] <- 1
                Y1[X1==1] <- 0
                X2       <- as.matrix(X2)
                Y2       <- X2
                Y2[X2==0] <- 1
                Y2[X2==1] <- 0
                return((t(X1) %*% X2 + t(Y1) %*% Y2)/normalization)
                }



RemoveSmallEigenvalues    <- function(M,min.val=0)  {
    EV      <- eigen(M)
    comps   <- which(EV$values<=min.val)
    M1      <- EV$vectors
    M1[,comps]<-0
    EV$vectors[comps]   <- 0
return(M1 %*% diag(EV$values) %*% t(M1))
}

# old :
ttest   <- function(snp.obj,trait.obj) {summary(lm(y ~ x, data=data.frame(x = matrix(as.numeric(snp.obj),ncol=1), y = as.numeric(trait.obj))))$coefficients[2,4]}

# after Stich (2008?)
IbdCorrection  <- function(K,T=0) {return(apply(K-matrix(T,nrow(K),ncol(K)),c(1,2),function(x) max(x,0))/(1-T))}

ComputeLargeCovarianceMatrix  <- function(mtr,file.name=NULL,block.size=50,correlation=TRUE) {
# INPUT :
# correlation : if TRUE, compute the correlation matrix; otherwise the covariance matrix
# OUTPUT :
# covariance matrix
    mtr  <- as.matrix(mtr)
    nind <- ncol(mtr)
    C    <- matrix(0,ncol(mtr),ncol(mtr))
    if (nind<=block.size) {
      C   <- cor(mtr)
      nbl <- 1
      blocks <- list(1:nind)
    } else {
      blocks          <- NULL
      nbl             <- ceiling(nind/block.size)
      if (nbl== nind/block.size) {
        for (i in 1:nbl) {blocks[[i]]   <- (1:nind)[(i-1)*block.size + 1:block.size]}
      } else {
        for (i in 1:(nbl-1)) {blocks[[i]]   <- (1:nind)[(i-1)*block.size + 1:block.size]}
        blocks[[nbl]]   <- (1:nind)[-(1:((nbl-1)*block.size))]
      }
    }

    for (i in 1:nbl) {
      if (correlation) {C[blocks[[i]],blocks[[i]]]  <- cor(mtr[,blocks[[i]]])} else {C[blocks[[i]],blocks[[i]]]  <- cov(mtr[,blocks[[i]]])}
        }
    if (nbl>1) {
      for (i in 2:nbl) {
        for (j in 1:i) {
          if (correlation) {C[blocks[[i]],blocks[[j]]]  <- cor(mtr[blocks[[i]]],mtr[,blocks[[j]]])} else {C[blocks[[i]],blocks[[j]]]  <- cov(mtr[blocks[[i]]],mtr[,blocks[[j]]])}
          C[blocks[[j]],blocks[[i]]]  <- t(C[blocks[[i]],blocks[[j]]])
        }
      }
    }
return(C)
}

#pst <- function(str) {paste(,sep="")}
# AddMeans(input.frame=GWAS.obj$pheno,col.select=3)[,ncol(GWAS.obj$pheno)+1]
AddMeans <- function(input.frame,col.select=1:ncol(input.frame),keep.original=TRUE) {
    # input.frame is assumed to have the same format as data1; it should contain a column named "genotype"
    #input.frame=data1;col.select=3:ncol(input.frame);keep.original=TRUE
    n.col         <- ncol(input.frame)
    n.col.select  <- length(col.select)
    input.frame   <- cbind(input.frame,input.frame[col.select])
    input.frame[,-(1:n.col)]  <- NA
    names(input.frame)[-(1:n.col)]  <- paste(names(input.frame)[col.select],".mean",sep="")
    #tapply(data1$AUC.NS, factor(data1$genotype,ordered=T),FUN=mean,na.rm =T)
    ind.names     <- unique(input.frame$genotype)
    first.occurence <- rep(0,length(ind.names))
    for (pl in 1:length(ind.names)) {first.occurence[pl]  <- min(which(input.frame$genotype==ind.names[pl]))}
    for (cl in 1:n.col.select) {
      new.col <- aggregate(input.frame[,col.select[cl]], by=list(input.frame$genotype),FUN=mean,na.rm =T)
      new.col <-     new.col[match(ind.names,new.col$Group.1),2]
      new.col[is.nan(new.col)]  <- NA
      input.frame[first.occurence,n.col+cl]   <- new.col
    }
    #if(!keep.original)
    return(input.frame)
}

Heritability <- function(data.vector,geno.vector)   {
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
# OUTPUT :
# heritability, genetic - and residual variance
    her.frame   <- data.frame(dat=data.vector,geno=geno.vector)
    her.frame   <- her.frame[!is.na(her.frame$dat),]
    her.frame$geno   <- factor(her.frame$geno)

    if (max(table(her.frame$geno))==1) {
      return(list(heritability=NA,gen.variance=NA,res.variance=NA))
    } else {
      av          <- anova(lm(dat~geno,data=her.frame))[[3]]
      return(list(heritability=av[1]/sum(av),gen.variance=av[1],res.variance=av[2]))
    }
}

#Heritability(data.vector=GWAS.obj$pheno[,tr.n],geno.vector=GWAS.obj$pheno$genotype)

Heritability2 <- function(data.vector,geno.vector)   {
# INCORRECT!!!!!!!!!!!!!!
#
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
# OUTPUT :
# heritability, genetic - and residual variance
    her.frame   <- data.frame(data=data.vector,geno=geno.vector)
    her.frame   <- her.frame[!is.na(her.frame$data),]
    means       <- aggregate(her.frame$data, by=list(her.frame$geno),FUN=mean,na.rm =T)[,2]
    vars        <- aggregate(her.frame$data, by=list(her.frame$geno),FUN=var,na.rm =T)[,2]
    return(list(heritability=var(means,na.rm =T)/(mean(vars,na.rm =T)+var(means,na.rm =T)),gen.variance=var(means,na.rm =T),res.variance=mean(vars,na.rm =T)))
}

Heritability3 <- function(data.vector,geno.vector,line.heritability=F,covariates.frame=data.frame())   {
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes. Must have the same length as geno.vector
#' @param line.heritability : if TRUE, the definition h2 = \sigma_g^2 / (\sigma_g^2 + \sigma_e^2 / r) is used, r being the number of replicates;
#                       otherwise (default) : h2 = \sigma_g^2 / (\sigma_g^2 + \sigma_e^2)
#' @param covariates.frame : a data.frame with additional covariates; must have the same number of rows as the length of data.vector and geno.vector
#                      Genotypes must also be in the same order as geno.vector
# OUTPUT :
# heritability, genetic - and residual variance
#
# data.vector=data[,tr.n] ; geno.vector= data[,1]; covariates.frame=data.frame(rep=data$Replicate)
# data.vector=data2[,j];geno.vector= data2[,1];covariates.frame=data.frame(rep=temp.rep.factor,x=data2$x_within_image,y=data2$y_within_image )
    stopifnot(length(geno.vector)==length(data.vector))
    her.frame   <- data.frame(dat=data.vector,geno=geno.vector)
    her.frame   <- her.frame[!is.na(her.frame$dat),]
    number.of.covariates <- 0

    if (nrow(covariates.frame) > 0) {
      stopifnot(nrow(covariates.frame)==length(data.vector))
      cov.names <- names(covariates.frame)
      number.of.covariates <- length(cov.names)
      covariates.frame <- as.data.frame(covariates.frame[!is.na(data.vector),])
      names(covariates.frame) <- cov.names
      her.frame <- cbind(her.frame,covariates.frame)
      her.frame <- her.frame[apply(covariates.frame,1,function(x){sum(is.na(x))})==0,]
      #cov.names <- names(her.frame)[-(1:2)]
    }

    her.frame$geno   <- factor(her.frame$geno)

    n.rep.vector <- as.integer(table(her.frame$geno))
    n.geno       <- length(n.rep.vector)

    # see Lynch and Walsh, Ch. 18, p. 559
    average.number.of.replicates <- ( sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector) ) / ( length(n.rep.vector) - 1 )

    if (max(n.rep.vector)==1) {
      return(list(heritability=NA,gen.variance=NA,res.variance=NA))
    } else {
      if (nrow(covariates.frame) > 0) {
        av          <- anova(lm(as.formula(paste('dat~geno+',paste(cov.names,collapse='+'))),data=her.frame))
      } else {
        av          <- anova(lm(dat~geno,data=her.frame))
      }
      gen.variance <- (av[[3]][1] - av[[3]][2+number.of.covariates]) / average.number.of.replicates
      if (gen.variance < 0) {gen.variance <- 0}
      if (line.heritability) {
        res.variance <- av[[3]][2+number.of.covariates] / average.number.of.replicates
      } else {
        res.variance <- av[[3]][2+number.of.covariates]
      }

    if (!line.heritability) {
      F.ratio        <- av[[3]][1] / av[[3]][2 + number.of.covariates]
      df.1           <- n.geno - 1
      #df.2           <- n.geno * (average.number.of.replicates - 1)
      df.2           <- av[[1]][2 + number.of.covariates]
      F.L            <- qf(0.025, df1=df.1, df2=df.2,lower.tail = TRUE)
      F.U            <- qf(0.975, df1=df.1, df2=df.2,lower.tail = TRUE)

      conf.int.left  <- (F.ratio / F.U - 1) / (F.ratio / F.U + average.number.of.replicates - 1)  # see Lynch and Walsh, Ch. 18, p. 563 (up to a factor 4)
      conf.int.right <- (F.ratio / F.L - 1) / (F.ratio / F.L + average.number.of.replicates - 1)
      if (conf.int.left < 0) {conf.int.left  <- 0}
      if (conf.int.right > 1) {conf.int.right  <- 1}
      if (conf.int.left > 1) {conf.int.left  <- 1}
      if (conf.int.right < 0) {conf.int.right  <- 0}

    } else { # to do...
      conf.int.left  <- NA
      conf.int.right <- NA
    }

    return(list(heritability=gen.variance/(gen.variance + res.variance),gen.variance=gen.variance,res.variance=res.variance,
                  line.heritability=line.heritability,average.number.of.replicates=average.number.of.replicates,
                  conf.int.left=conf.int.left,conf.int.right=conf.int.right))
    }
}

pin <- function (object, transform) {
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3)
  transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}

pin.adapted <- function (object, average.number.of.replicates=average.number.of.replicates, transform) {
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3)
  transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}


HeritabilityMixedModel <- function(data.vector,geno.vector,covariates=NULL,K,no.BLUPs=F,line.heritability=F,average.number.of.replicates=1)   {
# requires emma.R : emma.REMLE function
#
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
#' @param K           :  kinship matrix, of class matrix. Column-names should contain all occurring genotypes (?)
#' @param covariates  :  additional covariates (which should not include an intercept, which is included automatically)
#' @param line.heritability            :
#' @param average.number.of.replicates :

# OUTPUT :
# heritability, genetic - and residual variance
# data.vector=pheno.frame[,5]; geno.vector=pheno.frame$genotype; K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)]
    #data.vector=b$sim1;geno.vector=b$genotype

# data.vector=pheno.frame[,5]; geno.vector=pheno.frame$genotype; covariates=pheno.frame[,-(1:5)]; K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)]
#browser()
    #if (use.asreml) require(asreml)

#data.vector=GWAS.obj$pheno[,i];covariates=GWAS.obj$pheno[,c(2,5,6)];geno.vector=GWAS.obj$pheno$genotype;K=GWAS.obj$kinship;no.BLUPs=F;line.heritability=F;average.number.of.replicates=1
# data.vector=temp.pheno$sim;covariates=NULL;geno.vector=temp.pheno$genotype;K=K2[temp.pheno$genotype;temp.pheno$genotype];no.BLUPs=no.BLUPs;line.heritability=line.heritability;average.number.of.replicates=av.temp

    require(asreml)

    her.frame            <- data.frame(dat=data.vector,genotype=geno.vector)
    her.frame$genotype   <- as.character(her.frame$genotype)

    if (!is.null(covariates)) {
      covariates            <- as.data.frame(covariates)
      n.cov                 <- ncol(covariates)
      names(covariates)     <- paste('c',1:n.cov,sep='')
      row.names(covariates) <- row.names(her.frame)
      her.frame             <- cbind(her.frame,covariates)
    }

    her.frame     <- her.frame[!is.na(her.frame$dat),]
    n.rep.vector  <- as.integer(table(her.frame$genotype))

    # compute the average (or 'effective') number of replicates, see Lynch and Walsh, chapter 18
    average.number.of.replicates2 <- ( sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector) ) / ( length(n.rep.vector) - 1 )

    # If average.number.of.replicates is provided by the user and is larger than one, AND differs from average.number.of.replicates2 :
    # Then continue with the latter value
    if (average.number.of.replicates > 1 & average.number.of.replicates2 > 1) {
      average.number.of.replicates <- average.number.of.replicates2
    }

    if (line.heritability) {
      average.number.of.replicates <- average.number.of.replicates2
    }


    # re-scale the kinship matrix, for those genotypes for which phenotypic data are available
    K <- K[unique(her.frame$genotype),unique(her.frame$genotype)]
    K <- K / KinshipTransform(K)

    K.asreml <<- MakeKinshipAsreml(K)

    if (is.null(covariates)) {
      fixed.formula <<- as.formula('dat ~ 1')
    } else {
      cov.part <- paste(paste('c',1:n.cov,sep=''),collapse='+')
      fixed.formula <<- as.formula(paste('dat ~',cov.part))
    }

    random.formula <<- as.formula("~ giv(genotype,var=T)")
    #her.frame$genotype <- factor(her.frame$genotype)
    reml.obj <- asreml(fixed.formula,data=her.frame,random = random.formula,na.method.X="omit",ginverse = list(genotype=K.asreml))

    if (!no.BLUPs) {
      reml.pred <- predict(reml.obj, classify='genotype',data=her.frame)
      blup.frame <- reml.pred$predictions$pvals
      blup.frame$genotype <- as.character(blup.frame$genotype)
      row.names(blup.frame) <- blup.frame$genotype
    } else {
      blup.frame <- data.frame(genotype=unique(her.frame$genotype),predicted.value=rep(NA,length(unique(her.frame$genotype))))
      row.names(blup.frame) <- blup.frame$genotype
    }
    #head(blup.frame)

    var.comp <- summary(reml.obj)$varcomp$component

    inv.ai <- matrix(rep(0,4),ncol=2)
    inv.ai[lower.tri(inv.ai,diag=T)] <- reml.obj$ai
    inv.ai[1,2] <- inv.ai[2,1]

    if (line.heritability) {
        h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2] / average.number.of.replicates)
    } else {
        h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2] * average.number.of.replicates)
    }
    # in case of errors or strange output in predict.asreml, compute the BLUPs again
    if (!no.BLUPs) {
      if (max(abs(blup.frame$predicted.value),na.rm=T) > 100 * max(data.vector,na.rm=T)) {
        delta <- var.comp[1] / var.comp[2]
        K.temp <- K[her.frame$genotype,her.frame$genotype]
        V <- delta * K.temp + diag(ncol(K.temp))
        mu <- sum(solve(V) %*% as.matrix(her.frame$dat)) / sum(solve(V))
        Z.temp <- createZmatrixForFactor(f.vector=as.character(her.frame$genotype),factor.name="",ordered.factor=F)
        Z.temp <- Z.temp[,colnames(K)]
        blup.frame[colnames(Z.temp),]$predicted.value <- delta * K %*% t(Z.temp) %*% solve(V) %*% as.matrix(her.frame$dat-mu) + mu
      }
    }

    if (line.heritability) {
      st.error1 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ V1 / (V1 + V2/average.number.of.replicates))[2])
    } else {
      st.error1 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ V1 / (V1 + V2 * average.number.of.replicates))[2])
    }
    conf.int1 <- c(h2.estimate - qnorm(.975)*st.error1,h2.estimate + qnorm(.975)*st.error1)

    if (line.heritability) {
      st.error2 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ log(V2 / (V1*average.number.of.replicates)))[2])
    } else {
      st.error2 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ log(V2 * average.number.of.replicates / V1))[2])
    }

    if (line.heritability) {
      conf.int2 <- c(1 / (1 + exp(log(var.comp[2]/(var.comp[1]*average.number.of.replicates)) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[2]/(var.comp[1]*average.number.of.replicates)) - qnorm(.975)*st.error2)))
    } else {
      conf.int2 <- c(1 / (1 + exp(log(var.comp[2] * average.number.of.replicates /var.comp[1]) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[2]* average.number.of.replicates/var.comp[1]) - qnorm(.975)*st.error2)))
    }

    vg.conf.int.left   <- var.comp[1] - qnorm(0.975) * summary(reml.obj)$varcomp$std.error[1]
    vg.conf.int.right  <- var.comp[1] + qnorm(0.975) * summary(reml.obj)$varcomp$std.error[1]
    ve.conf.int.left   <- var.comp[2] - qnorm(0.975) * summary(reml.obj)$varcomp$std.error[2]
    ve.conf.int.right  <- var.comp[2] + qnorm(0.975) * summary(reml.obj)$varcomp$std.error[2]

return(list(vg=var.comp[1],ve=var.comp[2],h2=h2.estimate,conf.int1=conf.int1,conf.int2=conf.int2,inv.ai=inv.ai,
            vg.conf.int.left=vg.conf.int.left,vg.conf.int.right=vg.conf.int.right,ve.conf.int.left=ve.conf.int.left,ve.conf.int.right=ve.conf.int.right,
            blup.frame=blup.frame,line.heritability=line.heritability,average.number.of.replicates=average.number.of.replicates))
}

HeritabilityMixedModelWithFixedSigmaE2 <- function(data.vector,geno.vector,sigmaE2,covariates=NULL,K,no.BLUPs=F,line.heritability=F,average.number.of.replicates=1)   {
# requires emma.R : emma.REMLE function
#
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
#' @param K           :  kinship matrix, of class matrix. Column-names should contain all occurring genotypes (?)
# N.B. This is NOT the kinship matrix for the replicates
# OUTPUT :
# heritability, genetic - and residual variance
# data.vector=pheno.frame[,5]; geno.vector=pheno.frame$genotype; K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)]
    #data.vector=b$sim1;geno.vector=b$genotype

# data.vector=pheno.frame[,5]; geno.vector=pheno.frame$genotype; covariates=pheno.frame[,-(1:5)]; K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)]
#browser()
    #if (use.asreml) require(asreml)
# data.vector=temp.pheno$sim;geno.vector=temp.pheno$genotype;K=K2[temp.pheno$genotype,temp.pheno$genotype]; covariates=NULL; no.BLUPs=no.BLUPs; line.heritability=line.heritability;average.number.of.replicates=av.temp;sigmaE2=Heritability3(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)[3]
    require(asreml)
    # when using msm :
    #require(msm)
    #
    her.frame   <- data.frame(dat=data.vector,genotype=geno.vector)
    her.frame$genotype   <- as.character(her.frame$geno)

    if (!is.null(covariates)) {
      covariates <- as.data.frame(covariates)
      n.cov <- ncol(covariates)
      names(covariates) <- paste('c',1:n.cov,sep='')
      row.names(covariates) <- row.names(her.frame)
      her.frame <- cbind(her.frame,covariates)
    }

    her.frame   <- her.frame[!is.na(her.frame$dat),]
    #her.frame$geno   <- factor(her.frame$geno)
    #K <- K[levels(her.frame$geno),levels(her.frame$geno)]
    # corrected
    #n.rep.vector <- as.integer(table(her.frame$genotype))
    #n.rep.vector <- rep(n.rep.vector,times=n.rep.vector)

    #n.rep.vector <- table(her.frame$genotype)[her.frame$genotype]
    n.rep.vector <- as.integer(table(her.frame$genotype))

    #her.frame   <- data.frame(her.frame,wt=n.rep.vector)

    average.number.of.replicates2 <- ( sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector) ) / ( length(n.rep.vector) - 1 )
    if (average.number.of.replicates > 1 & average.number.of.replicates2 > 1) {average.number.of.replicates <- average.number.of.replicates2}

    K <- K[unique(her.frame$genotype),unique(her.frame$genotype)]
    #K <- K[her.frame$geno,her.frame$geno]
    K <- K / KinshipTransform(K)

    #if (use.asreml) {
    K.asreml <<- MakeKinshipAsreml(K)

    if (is.null(covariates)) {
      fixed.formula <<- as.formula('dat ~ 1')
    } else {
      cov.part <- paste(paste('c',1:n.cov,sep=''),collapse='+')
      fixed.formula <<- as.formula(paste('dat ~',cov.part))
    }
    random.formula <<- as.formula("~ giv(genotype,var=T)")
    #reml.obj <- asreml(sim1 ~ 1,data=b,random = random.formula,na.method.X="omit",ginverse = list(genotype=K.asreml))

    #reml.obj <- asreml(fixed.formula,data=her.frame,random = random.formula,na.method.X="omit",weights=wt,ginverse = list(genotype=K.asreml))
    reml.obj.init <- asreml(fixed.formula,data=her.frame,random = random.formula,na.method.X="omit",ginverse = list(genotype=K.asreml),rcov=~idv(units),start.values=T)
    iv            <- reml.obj.init$gammas.table
    iv$Value[3]  <- sigmaE2
    iv$Constraint[3]  <- 'F'
    #temp.init     <- sigmaE2
    #names(temp.init) <- 'F'
    #,rcov=~id(units,init=temp.init)
    reml.obj      <- asreml(fixed.formula,data=her.frame,random = random.formula,na.method.X="omit",ginverse = list(genotype=K.asreml),rcov=~idv(units),R.param=iv)

    if (!no.BLUPs) {
      reml.pred <- predict(reml.obj, classify='genotype',data=her.frame)
      blup.frame <- reml.pred$predictions$pvals
      blup.frame$genotype <- as.character(blup.frame$genotype)
      row.names(blup.frame) <- blup.frame$genotype
    } else {
      blup.frame <- data.frame(genotype=unique(her.frame$genotype),predicted.value=rep(NA,length(unique(her.frame$genotype))))
      row.names(blup.frame) <- blup.frame$genotype
    }
    #reml.obj <- asreml(as.formula('dat ~ 1'),data=her.frame,random = random.formula,na.method.X="omit",weights=wt,ginverse = list(genotype=solve(K)))
    var.comp <- summary(reml.obj)$varcomp$component

    # when using msm :
    inv.ai <- matrix(rep(0,4),ncol=2)
    inv.ai[lower.tri(inv.ai,diag=T)] <- reml.obj$ai
    inv.ai[1,2] <- inv.ai[2,1]

    #h2.estimate <- var.comp[1]/sum(var.comp)

    if (line.heritability) {
      #if (average.number.of.replicates > 1) {
        h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[3] / average.number.of.replicates)
      #} else {
      #  h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2])
      #}
    } else {
      #if (average.number.of.replicates > 1) {
      #  h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2])
      #} else {
        h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[3])# * average.number.of.replicates)
      #}
    }
    # in case of errors or strange output in predict.asreml, compute the BLUPs again
    if (!no.BLUPs) {
      if (max(abs(blup.frame$predicted.value),na.rm=T) > 100 * max(data.vector,na.rm=T)) {
        delta <- var.comp[1] / var.comp[3]
        K.temp <- K[her.frame$genotype,her.frame$genotype]
        V <- delta * K.temp + diag(ncol(K.temp))
        mu <- sum(solve(V) %*% as.matrix(her.frame$dat)) / sum(solve(V))
        Z.temp <- createZmatrixForFactor(f.vector=as.character(her.frame$genotype),factor.name="",ordered.factor=F)
        Z.temp <- Z.temp[,colnames(K)]
        #blup.frame[colnames(Z.temp),]$predicted.value <- delta * K.temp %*% t(Z.temp) %*% solve(V) %*% as.matrix(her.frame$dat-mu) + mu
        # corrected:
        blup.frame[colnames(Z.temp),]$predicted.value <- delta * K %*% t(Z.temp) %*% solve(V) %*% as.matrix(her.frame$dat-mu) + mu
      }
    }
    # when using msm :
    #st.error1 <- deltamethod( ~ x1 / (x1 + x2), mean=var.comp, cov=inv.ai, ses=TRUE)
    #cat(' qwedfsdfgbsdfgsdfg' )
    if (line.heritability) {
      #if (average.number.of.replicates > 1) {
      st.error1 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ V1 / (V1 + V3/average.number.of.replicates))[2])
      #} else {
      #  st.error1 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ V1 / (V1 + V2/average.number.of.replicates))[2])
      #}
    } else {
      #st.error1 <- as.numeric(pin(reml.obj,h2 ~ V1 / (V1 + V2))[2])
      st.error1 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ V1 / (V1 + V3 * average.number.of.replicates))[2])
    }
    #cat(' qwedfsdfgbsdfgsdfg5234523452345' )
    # old version :
    # st.error1 <- as.numeric(pin(reml.obj,h2 ~ V1 / (V1 + V2))[2])
    conf.int1 <- c(h2.estimate - qnorm(.975)*st.error1,h2.estimate + qnorm(.975)*st.error1)

    # when using msm :
    #st.error2 <- deltamethod( ~ log(x2 / x1), mean=var.comp, cov=inv.ai, ses=TRUE)
    # old version :
    # st.error2 <- as.numeric(pin(reml.obj,h2 ~ log(V2 / V1))[2])
    # conf.int2 <- c(1 / (1 + exp(log(var.comp[2]/var.comp[1]) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[2]/var.comp[1]) - qnorm(.975)*st.error2)))
    if (line.heritability) {
      st.error2 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ log(V3 / (V1*average.number.of.replicates)))[2])
    } else {
      st.error2 <- as.numeric(pin.adapted(reml.obj,average.number.of.replicates=average.number.of.replicates,h2 ~ log(V3 * average.number.of.replicates / V1))[2])
      #st.error2 <- as.numeric(pin(reml.obj,h2 ~ log(V2 / V1))[2])
    }
    if (line.heritability) {
      conf.int2 <- c(1 / (1 + exp(log(var.comp[3]/(var.comp[1]*average.number.of.replicates)) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[3]/(var.comp[1]*average.number.of.replicates)) - qnorm(.975)*st.error2)))
    } else {
      #conf.int2 <- c(1 / (1 + exp(log(var.comp[2]/var.comp[1]) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[2]/var.comp[1]) - qnorm(.975)*st.error2)))
      conf.int2 <- c(1 / (1 + exp(log(var.comp[3] * average.number.of.replicates /var.comp[1]) + qnorm(.975)*st.error2)),1 / (1 + exp(log(var.comp[3]* average.number.of.replicates/var.comp[1]) - qnorm(.975)*st.error2)))
    }
    #st.error2 <- deltamethod( ~ log(1 + x1 / x2), mean=var.comp, cov=inv.ai, ses=TRUE)
    #conf.int2 <- c(1 - exp(log(1+ var.comp[1]/var.comp[2]) - qnorm(.975)*st.error2),1 - exp(log(1+ var.comp[1]/var.comp[2]) + qnorm(.975)*st.error2))

    vg.conf.int.left   <- var.comp[1] - qnorm(0.975) * summary(reml.obj)$varcomp$std.error[1]
    vg.conf.int.right  <- var.comp[1] + qnorm(0.975) * summary(reml.obj)$varcomp$std.error[1]
    ve.conf.int.left   <- var.comp[2] - qnorm(0.975) * summary(reml.obj)$varcomp$std.error[2]
    ve.conf.int.right  <- var.comp[2] + qnorm(0.975) * summary(reml.obj)$varcomp$std.error[2]

return(list(vg=var.comp[1],ve=var.comp[3],h2=h2.estimate,conf.int1=conf.int1,conf.int2=conf.int2,inv.ai=inv.ai,
            vg.conf.int.left=vg.conf.int.left,vg.conf.int.right=vg.conf.int.right,ve.conf.int.left=ve.conf.int.left,ve.conf.int.right=ve.conf.int.right,
            blup.frame=blup.frame,line.heritability=line.heritability,average.number.of.replicates=average.number.of.replicates))
}


HeritabilityGE <- function(data.vector,geno.vector,env.vector)   {
# INPUT :
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
#' @param env.vector :  vector (character or factor) of environments
# OUTPUT :
# heritability over the environments
    her.frame   <- data.frame(pheno=data.vector,geno=geno.vector,env=env.vector)
    her.frame   <- her.frame[!is.na(her.frame$pheno),]
    her.frame$geno   <- factor(her.frame$geno)
    her.frame$env    <- factor(her.frame$env)
    av          <- anova(lm(pheno~env + geno + env:geno,data=her.frame))[[3]]
    return(av[2]/(av[2] + av[3]/nlevels(her.frame$env) +  av[1]/mean(table(her.frame$geno))))
}
#HeritabilityGE(data.vector=new.data$fw,geno.vector=new.data$genotype,env.vector=new.data$env)



GenomicControl <- function(LRT.stats) {
# OBJECTIVE : correction of LRT-stats based on the genomic inflation factor, as in Devlin and Roeder (1999?)
    inflation   <- median(LRT.stats)/0.456
    return(LRT.stats/inflation)
}

# OBJECTIVE :
# *
# *
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *

GenomicControlPvalues <- function(pvals,n.obs,n.cov=0) {
# OBJECTIVE : correction of p-values based on the genomic inflation factor, as in Devlin and Roeder (1999?)
# INPUT :
#' @param pvals : vector of p-values; may contain NA's
#' @param n.obs : number of individuals
#' @param n.cov : number of covariables
# OUTPUT :
# * the new, corrected p-values, with the same NA's as the input
# * genomic inflation factor
#
# assumes p-values from an F-test with df=1 and df2=n.obs-n.cov-2
    pvals[pvals==-1] <- 1  #  | is.na(pvals)
    new.pvals   <- pvals
    F.stats     <- qf(na.omit(pvals), df1=1, df2=n.obs-n.cov-2,lower.tail=F)
    inflation   <- median(F.stats,na.rm=T)/qf(0.5, df1=1, df2=n.obs-n.cov-2,lower.tail=F)
    F.stats     <- F.stats /inflation
    new.pvals[!is.na(pvals)]  <- pf(F.stats, df1=1, df2=n.obs-n.cov-2,lower.tail=F) # df2=400, previously ???
return(list(pvalues=new.pvals, inflation.factor=inflation))
}

GenomicControlPvaluesMTMM <- function(LOD.scores,n1,n2) {
# OBJECTIVE : correction of p-values based on the genomic inflation factor, as in Devlin and Roeder (1999?)
# INPUT :
#' @param LOD.scores : data-frame with 5 columns, containing the LOD-scores from MTMM
#        (LOD_Y1_univariate LOD_Y2_univariate    LOD_full LOD_trait_specific LOD_trait_common)
#
#' @param n1 : number of non-missing observations for trait 1
#' @param n2 : number of non-missing observations for trait 2

  new.LOD.scores <- LOD.scores
  df1 <- c(1,1,2,1,1)
  df2 <- c(n1-2,n2-2,n1+n2-4,n1+n2-4,n1+n2-3)

  for (j in 1:5) {
    markers.to.correct <- which(LOD.scores[,j]!=0 & !is.na(LOD.scores[,j]))
    pvals       <- 10^(-LOD.scores[markers.to.correct,j])
    F.stats     <- qf(pvals, df1=df1[j], df2=df2[j],lower.tail=F)
    inflation   <- median(F.stats,na.rm=T)/qf(0.5, df1=df1[j], df2=df2[j],lower.tail=F)
    F.stats     <- F.stats / inflation
    new.LOD.scores[markers.to.correct,j]  <- -log10(pf(F.stats, df1=df1[j], df2=df2[j],lower.tail=F))
  }
return(new.LOD.scores)
}


MakeSnpBoxplot    <- function(data.vector,marker.vector,file.name="") {
    data.vector   <- as.numeric(data.vector)
    marker.vector <- as.numeric(marker.vector)
    plot.data   <- data.frame(marker.value=marker.vector,trait.value=data.vector)
    plot.data   <- plot.data[!is.na(plot.data$trait.value),]
    if (file.name=="") {
        boxplot(trait.value ~ marker.value,data=plot.data,
                names=c(paste("n0=",as.character(sum(plot.data$marker.value==0)),sep=""),
                paste("n1=",as.character(sum(plot.data$marker.value==1)),sep="")))
    } else {
        jpeg(file.name,quality=100)
        boxplot(trait.value ~ marker.value,data=plot.data,
        names=c(paste("n0=",as.character(sum(plot.data$marker.value==0)),sep=""),
        paste("n1=",as.character(sum(plot.data$marker.value==1)),sep="")))
        dev.off()
    }
}
#make.snp.boxplot(data.vector=data1[,tr],snp.number=snp.sig[sg,1],file.name=paste(data.path,"boxplots/",trait,".snp.",as.character(snp.sig[sg,1]),Kmethod,suffix,".jpeg",sep=""))
# make.snp.boxplot(data.vector=data1[,tr],snp.number=snp.sig[37,1],file.name=paste(data.path,"plots/",trait,".snp.",as.character(snp.sig[37,1]),Kmethod,suffix,".jpeg",sep=""))

MakeSnpSummary    <- function(data.vector,trait.name="??",marker.vector,file.name="snp.summary.txt") {
# N.B. relies on the global objects snp and data1
    data.vector <- as.numeric(data.vector)
    marker.vector <- as.numeric(marker.vector)
    summary.data<- data.frame(genotype=data1$genotype,marker.value=marker.vector,trait.value=as.numeric(data.vector))
    summary.data<- summary.data[!is.na(summary.data$trait.value),]
    cat("Trait: ",trait.name,"\n",file=file.name)
    cat("snp number: ",snp.number,"\n","\n",file=file.name,append=TRUE)
    cat("Number of accessions with at least one nonmissing observation: ",length(unique(as.character(summary.data$genotype))),"\n",file=file.name,append=TRUE)
    cat("Total number of nonmissing observations: ",nrow(summary.data),"\n","\n",file=file.name,append=TRUE)
    cat("Proportion Columbia allele at accession level: ",mean(as.numeric(snp[snp.number,unique(as.character(summary.data$genotype))])),"\n",file=file.name,append=TRUE)
    cat("Proportion Columbia allele at plant level: ",mean(summary.data$marker.value),"\n","\n",file=file.name,append=TRUE)
    #
    summary.data0           <- summary.data[summary.data$marker.value==0,]
    summary.data1           <- summary.data[summary.data$marker.value==1,]
    #
    suppressWarnings(kt0 <- ks.test(x=summary.data0$trait.value,y="pnorm"))
    suppressWarnings(kt1 <- ks.test(x=summary.data1$trait.value,y="pnorm"))
    #
    cat("p-value of Kolmogorov-Smirnov test for normality, for the plants with the 0 allele: ",kt0$p.value,"\n",file=file.name,append=TRUE)
    cat("p-value of Kolmogorov-Smirnov test for normality, for the plants with the 1 allele: ",kt1$p.value,"\n","\n",file=file.name,append=TRUE)
    #
    outliers0   <- outlier.test(summary.data0$trait.value,2)
    outliers1   <- outlier.test(summary.data1$trait.value,2)
    #
    if (length(outliers0$outliers)>0) {
            summary.data0           <- summary.data0[outliers0$outliers,]
            #rare.genotypes          <- unique(as.character(summary.data0$genotype))
            # aggregate(summary.data0$trait.value, by=list(summary.data0$genotype),FUN=length)
            cat("Outliers in genotypes with allele ",0," :","\n","\n",file=file.name,append=TRUE)
            suppressWarnings(write.table(summary.data0[,c(1,3)],quote=F,file=file.name,append=TRUE,row.names=F))
    }
    if (length(outliers1$outliers)>0) {
            summary.data1           <- summary.data1[outliers1$outliers,]
            cat("\n","Outliers in genotypes with allele ",1," :","\n","\n",file=file.name,append=TRUE)
            suppressWarnings(write.table(summary.data1[,c(1,3)],quote=F,file=file.name,append=TRUE,row.names=F))
    }
# to do : look up standars methods for outlier-detection
# insert a data.frame with columns : genotype, total.number.of.nonmissing, number.of.outliers
# also : outliers of the residuals ?
# also : output on screen as option
}


MakeHaploBoxplot    <- function(data.vector,hap.vector,file.name="",haplo.names=NULL,y.lab=NULL) {
    plot.data   <- data.frame(trait.value=as.numeric(data.vector),hap.vector=as.factor(hap.vector))
    if (file.name=="") {
      boxplot(trait.value ~ hap.vector,data=plot.data,names=haplo.names,ylab=y.lab) # as.character(tabulate(hap.vector))
    } else {
      jpeg(file.name,quality=100)
      boxplot(trait.value ~ hap.vector,data=plot.data,names=haplo.names,ylab=y.lab) # as.character(tabulate(hap.vector))
      dev.off()
    }
}

# For an extended version including a training and validation set, see the function define.haplotypes in the SHARE3 script
FindHaplotypes <- function(marker.frame,snp.set=1:nrow(marker.frame),
                           ind.set=1:ncol(marker.frame),MAF=0.01,classification.constant=1) {
    # marker.frame is assumed to be a data-frame or matrix with the snps in the rows and individuals in the columns
    # The species is assumed to be homozygote (this function is an adapted version of the more general one in the haplo.stats package,
    # which also considers heterozygozity)
    phased.snp.data <- t(marker.frame[snp.set,ind.set])
    nr          <- nrow(phased.snp.data)
    hapSeq      <- apply(phased.snp.data, 1, function(x) {paste(x, sep = "", collapse = "-")})
    uniHap      <- unique(hapSeq)
    nHap        <- length(uniHap)
    hapNum      <- as.character(1:nHap) # N.B. This line replaces the following line from the original code (from the haplo.stats package) # ok now???
    #hapNum      <- unlist(sapply(1:nHap, function(x) {paste(paste(rep("0", ceiling(log10(nHap)) - nchar(as.character(x))),collapse = ""), x, sep = "",collapse = "")}))
    names(uniHap) <- paste("hap.", hapNum, sep = "")
    hapPool     <- strsplit(uniHap, "-")
    hapPool     <- data.frame(t(data.frame(lapply(hapPool, function(x) {as.numeric(x)}))))
    colnames(hapPool) <- colnames(phased.snp.data)
    hapCount    <- sapply(uniHap, function(x) {sum(x == hapSeq)})
    hapFreq     <- hapCount/sum(hapCount)
    noNameUniHap <- uniHap
    names(noNameUniHap) <- NULL
    hapIndex    <- sapply(hapSeq, function(x) {which(x == noNameUniHap)})
    hapIndex2   <- hapIndex
    # Now define respectively haplotypes which
    # 3) the subset of tr.types that has minimal (MAF) frequency (within tr.types, not the whole sample)
    # 4) the subset of tr.types that is below this frequency
    common.types<- sort((1:nHap)[(tabulate(hapIndex,nbins=nHap)/nr)>MAF])
    rare.types  <- sort((1:nHap)[(tabulate(hapIndex,nbins=nHap)/nr)<=MAF])
    #
    if (length(common.types)==0) {
      common.types <- rare.types
      rare.types <- numeric(0)
      no.common.type <- TRUE
    } else {
      no.common.type <- FALSE
    }
    hapIndex2[hapIndex %in% rare.types] <- 0
    if ((length(rare.types))>0) {
        new.hapPool     <- hapPool[common.types,]
        new.nHap        <- nrow(new.hapPool)
        new.hapIndex    <- matrix(0,nr,length(common.types))
        common.individuals<- which(hapIndex %in% common.types)
        n.c             <- length(common.individuals)# as.numeric(as.factor(... : this is to rename the common types
        new.hapIndex[matrix(c(common.individuals, as.numeric(as.factor(hapIndex[common.individuals]))), ncol=2)]   <- 1
        for (h in rare.types) {
            qw          <- hapPool[h,]
            differences <- apply(hapPool[common.types,],1,FUN=function(x) {sum(x!=qw)})
            weights     <- exp(-classification.constant*differences)
            weights     <- weights/sum(weights)
            new.hapIndex[which(hapIndex==h),]<- matrix(rep(weights,sum((hapIndex==h))),byrow=TRUE,ncol=length(weights))
        }
    hapPool     <- new.hapPool
    hapIndex     <- new.hapIndex
    nHap        <- new.nHap
    } else {
        new.hapIndex    <- matrix(0,nr,length(common.types))
        common.individuals<- which(hapIndex %in% common.types)
        n.c             <- length(common.individuals)# as.numeric(as.factor(... : this is to rename the common types
        new.hapIndex[matrix(c(common.individuals, as.numeric(as.factor(hapIndex[common.individuals]))), ncol=2)]   <- 1
        hapIndex     <- new.hapIndex
    }
return(list(hapIndex=hapIndex,nHap=nHap,hapPool=hapPool,hapCodes=hapIndex2,no.common.type=no.common.type))
}

# then :
#- asreml # nb columbia genotypes !!
#- call this using snp.sig, + genes
# license : system command !
# boxplots
# output in files
# document denife.haplotypes; classification, MAF etc
# interpretation haplo-tested
# correction haplo boxplot (when MAF>0)

############################################################################################################################


fdp <- function(pvals,ev=1,gamma=0.05) {
# computes false discovery proportion under the full (strong) null
p<- length(pvals)
if (gamma==0) {
    lod.thr= -log(ev/p,base=10)
    R.thr= sum(pvals<= ev/p)
} else {
    evtable = data.frame(pvalues=sort(pvals),test=rep(0,p))
    evtable$test  = as.numeric((p*evtable$pvalues) /(1:p)<= gamma)
    if (sum(evtable$test)==p) {
      R.thr=p
      lod.thr=-log(evtable$pvalues[p],base=10)
    } else {
      ind = min(which(evtable$test==0))-1
      if (ind==0) {R.thr=0; lod.thr=Inf} else {R.thr=ind; lod.thr=-log(evtable$pvalues[ind],base=10)}
    }
}
list(lod.thr,R.thr)
}

# estimate.interaction
#geno.vector=GWAS.obj$pheno$genotype;marker1=as.numeric(GWAS.obj$markers[interactions$snp1[ep],GWAS.obj$pheno$genotype]);marker2=as.numeric(GWAS.obj$markers[interactions$snp2[ep],GWAS.obj$pheno$genotype]);trait.values=GWAS.obj$pheno[,tr];cov.values=cov.frame.interaction;kinship.asreml.object=GWAS.obj$kinship.asreml
#
#
# assumes 0-1 markers
# TO DO : HETEROZYGOTES
#

EstimateInteraction  <- function(geno.vector,marker1,marker2,trait.values,kinship.asreml.object,cov.values=data.frame(NULL)) {
  #
  require(asreml)
  # N.B. the input geno.vector is (supposed to be) character; for asreml, ordered factor is requires
  interaction.frame   <- data.frame(genotype=ordered(geno.vector,levels=unique(geno.vector)),phenotype=trait.values,marker1=marker1,marker2=marker2,epistasis=marker1*marker2)
  if (ncol(cov.values)>1) {interaction.frame   <- cbind(interaction.frame,cov.values)}
  if (ncol(cov.values)==1) {
          reml.formula.inter        <- as.formula("phenotype ~ marker1 + marker2 + epistasis")
      } else {
          reml.formula.inter        <- as.formula(paste("phenotype ~ marker1 + marker2 + epistasis +",paste(names(cov.values),collapse="+")))
      }
  reml.obj<- asreml(fixed= reml.formula.inter,data=interaction.frame, random = ~ giv(genotype,var=T),na.method.X="omit",
                    ginverse = list(genotype=GWAS.obj$kinship.asreml),G.param=iv,R.param=iv,fixgammas=TRUE)
  # N.B. in the preceding line, we should have kinship.asreml.object instead of GWAS.obj$kinship.asreml; but the former gives an error. WHY ??
  # N.B another "global" object is used: iv, defined in emmax_win.r
  #  this is done to pass on the values of the variance components as they are in the scan for the main effects
  effect  <-  reml.obj$coefficients$fixed[["epistasis"]]
  pvalue  <- (wald(reml.obj))[[4]][4]
return(list(pvalue,effect))
}
#estimate.interaction(geno.vector=GWAS.obj$pheno$genotype,marker1=as.numeric(GWAS.obj$markers[interactions$snp1[ep],GWAS.obj$pheno$genotype]),marker2=as.numeric(GWAS.obj$markers[interactions$snp2[ep],GWAS.obj$pheno$genotype]),trait.values=GWAS.obj$pheno[,tr],kinship.asreml.object=GWAS.obj$kinship.asreml)

OutlierTest    <- function(x,number.of.sds=3) {
        st.dev  <- sd(x,na.rm=TRUE)
        avg     <- mean(x,na.rm=TRUE)
        outliers<- which(x < (avg - number.of.sds*st.dev) | x > (avg + number.of.sds*st.dev))
return(list(outlier.values=x[outliers],outliers=outliers))
}
#x       <- rnorm(1000)
#qw      <- outlier.test(x,6)


"stabilityselection" <- function(x,y,nbootstrap=100,nsteps=20,alpha=0.2,plotme=FALSE)
{
	# Stability selection in the spirit of Meinshausen&Buhlman
	# JP Vert, 14/9/2010

	# x is the n*p design matrix, y the n*1 variable to predict
	# x should be normalized to unit variance per column before calling this function (if you want to)
	# the result is a score (length p) for each feature, the probability that each feature is selected
  # during the first nsteps steps of the Lasso path when half of the samples are used and the features
  # are reweigthed by a random weight uniformaly sampled in [alpha,1].
  # This probability is estimated by nbootstrap bootstrap samples
	require(lars)
	dimx <- dim(x)
	n <- dimx[1]
	p <- dimx[2]
	halfsize <- as.integer(n/2)
	freq <- matrix(0,nsteps+1,p)

	for (i in seq(nbootstrap)) {

		# Randomly reweight each variable
		xs <- t(t(x)*runif(p,alpha,1))

		# Ramdomly split the sample in two sets
		perm <- sample(dimx[1])
		i1 <- perm[1:halfsize]
		i2 <- perm[(halfsize+1):n]

		# run the randomized lasso on each sample and check which variables are selected
		r <- lars(xs[i1,],y[i1],max.steps=nsteps,normalize=FALSE,use.Gram=FALSE)
		freq <- freq + abs(sign(coef.lars(r)))
		r <- lars(xs[i2,],y[i2],max.steps=nsteps,normalize=FALSE,use.Gram=FALSE)
		freq <- freq + abs(sign(coef.lars(r)))
		}

	# normalize frequency in [0,1]
	freq <- freq/(2*nbootstrap)

	if (plotme) {
		matplot(freq,type='l',xlab="LARS iteration",ylab="Frequency")
	}

	# the final stability score is the maximum frequency over the steps
	result <- apply(freq,2,max)
}

# END of code by JP Vert


DefineBlocks  <- function(indices,block.size=1000) {
    nsnp            <- length(indices)
    if (nsnp<=block.size) {
      blocks  <- list(indices)
    } else {
      blocks          <- NULL
      nbl             <- ceiling(nsnp/block.size)
      if (nbl== nsnp/block.size) {
          for (i in 1:nbl) {blocks[[i]]   <- indices[(i-1)*block.size + 1:block.size]}
      } else {
          for (i in 1:(nbl-1)) {blocks[[i]]   <- indices[(i-1)*block.size + 1:block.size]}
          blocks[[nbl]]   <- indices[-(1:((nbl-1)*block.size))]
      }
    }
return(blocks)
}


MakeKinshipAsreml  <- function(K,genotype.names=colnames(K)) {
# K = kinship matrix
    require(MASS)
    n                   <- nrow(K)
    vec1                <- rep(1:n,1:n)
    row.number.matrix   <- matrix(rep(1:n,n),ncol=n)
    vec2                <- row.number.matrix[upper.tri(row.number.matrix,diag=T)]
    matrix.indices      <- matrix(c(vec1,vec2),ncol=2)
    #
    Ainv        <- ginv(K)
    AINV        <- data.frame(matrix(0,round(.5*n*(n+1)),3))
    names(AINV) <- c("Row","Column","Ainverse")
    AINV[,1]    <- vec1
    AINV[,2]    <- vec2
    AINV[,3]    <- Ainv[matrix.indices]
    attr(AINV,"rowNames") <-  genotype.names
return(AINV )
}

MakeVarcompFile <- function(file.name,var.comp.values=c(1,0.1)) {
    var.comp.values      <- data.frame(var.comp=as.numeric(unlist(var.comp.values)))
    row.names(var.comp.values)<- c("sigma2_g","sigma2_e")
    write.table(var.comp.values,file=file.name,quote=F,row.names = T,col.names=F,sep=",")
}

MakePhenoFile <- function(pheno.object,file.name,col.number=2) {
    pheno.frame <- data.frame(pheno.object$genotype,pheno.object[,col.number])
    names(pheno.frame)  <- c("genotype",names(pheno.object)[col.number])
    write.csv(pheno.frame,quote=F,file=file.name)
}

MakeCovariateFile <- function(pheno.dataframe,cov.cols,file.name="") {
        cov.frame               <- data.frame(mu=rep(1,nrow(pheno.dataframe)),row.names=row.names(pheno.dataframe))
        if (sum(cov.cols)!=0) {
            new.names <- c(names(cov.frame),names(pheno.dataframe)[cov.cols])
            cov.frame <- cbind(cov.frame,pheno.dataframe[,cov.cols])
            names(cov.frame)  <- new.names
            write.csv(cov.frame,quote=F,row.names=F,file=file.name)
        }
return(cov.frame)
}

# to do: remove the temporary files after the bin-file has been created
MakeBin  <- function(kinship,plant.names,pheno,csv.file.name,bin.file.name) {
# old argument: markers : removed
    varcompfile <- "temp.varcomp.csv"
    pheno <- data.frame(genotype=plant.names,pheno=rnorm(length(plant.names)))
    vcov.matrix <- MakeScanGlsKinship(1,1,kinship,plant.names=plant.names,pheno,tr.n=2)
    write.table(vcov.matrix,file=varcompfile,sep=",",quote=F,col.names=F,row.names=F)
    #MakeVarcompFile(file.name=varcompfile)
    temp.pheno.file <- "temp.pheno.csv"
    MakePhenoFile(pheno.object=pheno,col.number=2,file.name=temp.pheno.file)
    output.file <- "temp.output.csv"
    command.string      <- paste("scan_GLS",csv.file.name,temp.pheno.file,varcompfile,output.file,"-writebin",bin.file.name)
    cat(command.string,"\n")
    system(command.string, intern = TRUE, ignore.stderr = FALSE,wait = TRUE, input = NULL)
}

# to do: remove the temporary files after the bin-file has been created
MakeBinOld0  <- function(kinship,plant.names,pheno,csv.file.name,bin.file.name) {
# old argument: markers : removed
    varcompfile <- "temp.varcomp.csv"
    pheno[,2] <- rnorm(nrow(pheno))
    vcov.matrix <- MakeScanGlsKinship(1,1,kinship,plant.names=plant.names,pheno,tr.n=2)
    write.table(vcov.matrix,file=varcompfile,sep=",",quote=F,col.names=F,row.names=F)
    #MakeVarcompFile(file.name=varcompfile)
    temp.pheno.file <- "temp.pheno.csv"
    MakePhenoFile(pheno.object=pheno,col.number=2,file.name=temp.pheno.file)
    output.file <- "temp.output.csv"
    command.string      <- paste("scan_GLS",csv.file.name,temp.pheno.file,varcompfile,output.file,"-writebin",bin.file.name)
    cat(command.string,"\n")
    system(command.string, intern = TRUE, ignore.stderr = FALSE,wait = TRUE, input = NULL)
}


#Old version", for scan_GLS prior to version 1.17 :
# to do: remove the temporary files after the bin-file has been created
MakeBinOld  <- function(markers,kinship.file,csv.file.name,bin.file.name) {
    varcompfile <- "temp.varcomp.csv"
    MakeVarcompFile(file.name=varcompfile)
    temp.pheno.file <- "temp.pheno.csv"
    MakePhenoFile(pheno.object=data.frame(genotype=rep(names(markers),each=2),temp.trait=rnorm(2*ncol(markers))),col.number=2,file.name=temp.pheno.file)
    output.file <- "temp.output.csv"
    command.string      <- paste("scan_GLS",csv.file.name,temp.pheno.file,kinship.file,varcompfile,output.file,"-writebin",bin.file.name)
    system(command.string, intern = TRUE, ignore.stderr = FALSE,wait = TRUE, input = NULL)
}

MakeCsv  <- function(markers,plant.names,file.name) {
    blocks <- DefineBlocks(1:nrow(markers),10000)
    cat("",plant.names,"\n",file=file.name,sep=",")
    for (i in 1:length(blocks)) {write.table(markers[blocks[[i]],],file=file.name,quote=F,append=T,col.names =F,sep=",")}
}

###########################
ReadMarkerAndGeneData  <- function(gffFileName="TAIR10_GFF3_genes.gff",geno.data="atwell_data.RData",
                                   snp.file.name="call_method_32_without_header.csv",
                                   header.name="call_method_32_header.csv",n.snp=216130,
                                   snp.start.col=3,ref.nr=1,maf=0.10) {
### OBJECTIVE :
#
### INPUT :
# gffFileName   : gff(3) file obtained from TAIR
# snp.file.name : the file containing all snp-data, WITHOUT header
#                  Should contain letters A,C,G,T and "-" for missing values
# geno.data     : name of the RData file that is to be created
# header.name   : file containing the accession names (see the default for an example)
#                 SHOULD JUST CONTAIN THE ACCESSION NAMES, ON ONE LINE, SEPARATED BY COMMAS!
# n.snp         : number of markers in the file snp.file.name
#                 since snp.file.name does not have a header (i.e. column names)
#                 n.snp should be the number of lines in this file.
#                 n.snp may be overestimated by several lines, but not underestimated
# snp.start.col : column number where the marker-data start. SHOULD NOW BE 3 !!
# ref.nr        : column number of the reference genotype, counted from snp.start.col
#                  i.e. if the first column of markers is the reference genotype, then ref.nr=1
# maf : minor allele frequency to which the markers in marker.object WILL BE restricted
#
### OUTPUT :
# an R-image with the name r.image.name, containing a data-frame snp, the numbers n (of genotypes) and N (of snps)
# in the data-frame, and the file plant.names.csv, containing all genotype names
# The data frame snp has N rows and n+2 columns. The first and second  column  contain the chromosome number and base pair
# position of each marker.
### comments

gc()
#header          <- read.table(file=header.name,sep=",")
#plant.names     <- as.character(unlist(header[1,snp.start.col:ncol(header)]))
plant.names     <- as.character(read.table(file=header.name,sep=",",colClasses="character"))
n               <- length(plant.names)                     # number of genotypes

list.of.blocks  <- DefineBlocks(1:n.snp,block.size=40000)
n.block         <- length(list.of.blocks)
con             <- file(description=snp.file.name, open = "r")
cat("",file=paste("new_",snp.file.name,sep=""))

for (bl in 1:n.block) {
  snp             <- read.table(con,sep=",", nrows = length(list.of.blocks[[bl]])) # ,colClasses=c(rep("character",n+2))
  number.of.missing <- apply(snp,1,FUN=function(x){sum(x=="-")+sum(is.na(x))})
  # remove the markers with at least one missing obervation
  snp             <- snp[number.of.missing==0,]
  if (maf>0) {
    minor.allele.fr <- apply(snp[,snp.start.col:ncol(snp)],1,FUN=function(x){sum(sort(table(unlist(x)),decreasing=T)[-1])})/n
    snp             <- snp[minor.allele.fr>=maf,]
  }
  snp[,snp.start.col:ncol(snp)] <- as.data.frame(lapply(snp[,snp.start.col:ncol(snp)],FUN=ConvertACGTfactorLevelsTo1234))
  # for all columns of the snp-data, except the column corresponding to the reference genotype,
  # we change the snp-values : they will be zero if a genotype is the same as the reference genotype
  # on that locus, and larger than zero otherwise
  snp[,-c(1:(snp.start.col-1),snp.start.col+ref.nr-1)] <- as.data.frame(lapply(snp[,-c(1:(snp.start.col-1),snp.start.col+ref.nr-1)],FUN=function(x){abs(x-snp[,snp.start.col+ref.nr-1])}))
  # We now put the snp-values of the reference genotype to zero
  snp[,snp.start.col+ref.nr-1] <- 0
  snp[,snp.start.col:ncol(snp)] <- as.data.frame(lapply(snp[,snp.start.col:ncol(snp)],FUN=function(x){x[x>1]<-1;return(1-x)}))
  # Keep the snp-data abd first two columns, which are supposed to contain the chromosome and position data
  # Throw away all other columns
  snp <- snp[,c(1:2,snp.start.col:ncol(snp))]
  write.table(snp,sep=",",col.names=F,row.names=F,quote=F,file=paste("new_",snp.file.name,sep=""),append=TRUE)
}
close(con)
# snp             <- read.table(file="large_test_file.csv",colClasses=rep("integer",10),sep=",")

rm(snp)
gc()

snp             <- read.table(file=paste("new_",snp.file.name,sep=""),sep=",",colClasses=rep("integer",10))
names(snp)      <- c("chromosome","position",plant.names)
row.names(snp)  <- paste("m",as.character(1:nrow(snp)),sep="")
N               <- nrow(snp)
#


gff <- gffRead(gffFileName)
gff <- gff[gff$feature=="gene",]
gff <- gff[(gff$seqname!="ChrC") & (gff$seqname!="ChrM"),]
gff$seqname     <- as.integer(as.factor(gff$seqname))
gff$Name <- getAttributeField(gff$attributes, "Name")

gc()

nchr            <- length(unique(snp[,1]))
chr.lengths     <- as.integer(table(snp[,1]))
chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])

gene.lengths     <- as.integer(table(gff$seqname))
gene.pos         <- c(0,cumsum(gene.lengths)[1:(nchr-1)])

for (i in 1:nchr) {
  gff$end[gene.pos[i] + 1:(gene.lengths[i]-1)]  <-  gff$start[gene.pos[i] + 2:(gene.lengths[i])] - 1
}

gene.info       <- rep(NA,nrow(snp))

for (i in 1:nrow(gff)) {
  gene.info[(snp$chromosome==gff$seqname[i]) & (snp$position %in% ((gff$start[i]):(gff$end[i])))]  <- gff$Name[i]
  #cat(i,"\n")
}

n.g         <- length(unique(gene.info[!is.na(gene.info)]))
genes       <- data.frame(gene.name=as.character(unique(gene.info[!is.na(gene.info)])),first.marker=rep(0,n.g),last.marker=rep(0,n.g),gene.length=rep(0,n.g))
for (i in 1:n.g) {gene.i  <- which(gene.info==genes$gene.name[i])
                  genes$first.marker[i] <- min(gene.i)
                  genes$last.marker[i]  <- max(gene.i)
                  }
genes$gene.length <- genes$last.marker-genes$first.marker+1
genes$gene.name   <-  as.character(genes$gene.name)

save(genes,n.g,gene.info,snp,n,N,plant.names,maf,file=geno.data)

} # end of ReadMarkerAndGeneData function

############################

# new version
ReadMarkerAndGeneData2  <- function(gffFileName="TAIR10_GFF3_genes.gff",geno.data="atwell_data.RData",
                                   snp.file.name="call_method_32_without_header.csv",
                                   header.name="call_method_32_header.csv",n.snp=216130,
                                   snp.start.col=3,ref.nr=1,maf=0.10,max.shift=3000) {
### OBJECTIVE :
#
### INPUT :
# gffFileName   : gff(3) file obtained from TAIR
# snp.file.name : the file containing all snp-data, WITHOUT header
#                  Should contain letters A,C,G,T and "-" for missing values
# geno.data     : name of the RData file that is to be created
# header.name   : file containing the accession names (see the default for an example)
#                 SHOULD JUST CONTAIN THE ACCESSION NAMES, ON ONE LINE, SEPARATED BY COMMAS!
# n.snp         : number of markers in the file snp.file.name
#                 since snp.file.name does not have a header (i.e. column names)
#                 n.snp should be the number of lines in this file.
#                 n.snp may be overestimated by several lines, but not underestimated
# snp.start.col : column number where the marker-data start. SHOULD NOW BE 3 !!
# ref.nr        : column number of the reference genotype, counted from snp.start.col
#                  i.e. if the first column of markers is the reference genotype, then ref.nr=1
# maf : minor allele frequency to which the markers in marker.object WILL BE restricted
# max.shift : maximum number of base pairs the end of a gene is extended to the right (in case the strand is -)
#                                          the start of a gene is extended to the left (in case the strand is +)
### OUTPUT :
# an R-image with the name r.image.name, containing a data-frame snp, the numbers n (of genotypes) and N (of snps)
# in the data-frame, and the file plant.names.csv, containing all genotype names
# The data frame snp has N rows and n+2 columns. The first and second  column  contain the chromosome number and base pair
# position of each marker.
### comments

gc()
#header          <- read.table(file=header.name,sep=",")
#plant.names     <- as.character(unlist(header[1,snp.start.col:ncol(header)]))
plant.names     <- as.character(read.table(file=header.name,sep=",",colClasses="character"))
n               <- length(plant.names)                     # number of genotypes

list.of.blocks  <- DefineBlocks(1:n.snp,block.size=40000)
n.block         <- length(list.of.blocks)
con             <- file(description=snp.file.name, open = "r")
cat("",file=paste("new_",snp.file.name,sep=""))

for (bl in 1:n.block) {
  snp             <- read.table(con,sep=",", nrows = length(list.of.blocks[[bl]])) # ,colClasses=c(rep("character",n+2))
  number.of.missing <- apply(snp,1,FUN=function(x){sum(x=="-")+sum(is.na(x))})
  # remove the markers with at least one missing obervation
  snp             <- snp[number.of.missing==0,]
  if (maf>0) {
    minor.allele.fr <- apply(snp[,snp.start.col:ncol(snp)],1,FUN=function(x){sum(sort(table(unlist(x)),decreasing=T)[-1])})/n
    snp             <- snp[minor.allele.fr>=maf,]
  }
  snp[,snp.start.col:ncol(snp)] <- as.data.frame(lapply(snp[,snp.start.col:ncol(snp)],FUN=ConvertACGTfactorLevelsTo1234))
  # for all columns of the snp-data, except the column corresponding to the reference genotype,
  # we change the snp-values : they will be zero if a genotype is the same as the reference genotype
  # on that locus, and larger than zero otherwise
  snp[,-c(1:(snp.start.col-1),snp.start.col+ref.nr-1)] <- as.data.frame(lapply(snp[,-c(1:(snp.start.col-1),snp.start.col+ref.nr-1)],FUN=function(x){abs(x-snp[,snp.start.col+ref.nr-1])}))
  # We now put the snp-values of the reference genotype to zero
  snp[,snp.start.col+ref.nr-1] <- 0
  snp[,snp.start.col:ncol(snp)] <- as.data.frame(lapply(snp[,snp.start.col:ncol(snp)],FUN=function(x){x[x>1]<-1;return(1-x)}))
  # Keep the snp-data abd first two columns, which are supposed to contain the chromosome and position data
  # Throw away all other columns
  snp <- snp[,c(1:2,snp.start.col:ncol(snp))]
  write.table(snp,sep=",",col.names=F,row.names=F,quote=F,file=paste("new_",snp.file.name,sep=""),append=TRUE)
}
close(con)
# snp             <- read.table(file="large_test_file.csv",colClasses=rep("integer",10),sep=",")

rm(snp)
gc()

snp             <- read.table(file=paste("new_",snp.file.name,sep=""),sep=",",colClasses=rep("integer",10))
names(snp)      <- c("chromosome","position",plant.names)
row.names(snp)  <- paste("m",as.character(1:nrow(snp)),sep="")
N               <- nrow(snp)
#
nchr            <- length(unique(snp[,1]))
chr.lengths     <- as.integer(table(snp[,1]))
chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])

### read gff file
gff <- gffRead(gffFileName)
gff <- gff[gff$feature %in% c("gene","transposable_element_gene","pseudogene"),]
gff <- gff[(gff$seqname!="ChrC") & (gff$seqname!="ChrM"),]
gff$seqname     <- as.integer(as.factor(gff$seqname))
gff$Name <- getAttributeField(gff$attributes, "Name")

### remove genes that are completely covered (overlapped) by another gene

geneRegion <- function(currentGeneNumber,End,Start=1,size=5) {max(Start,currentGeneNumber-size):min(End,currentGeneNumber+size)}
n.gene     <- nrow(gff)
gene.test  <- rep(0,n.gene)
for (i in 1:n.gene) {gene.test[i] <- sum(apply(gff[geneRegion(currentGeneNumber=i,End=n.gene),4:5],1,FUN=function(x){sum(gff[i,4:5] %in% (x[1]):(x[2]))})==2)}
gff        <- gff[which(gene.test==1),]

### sort gff (with starting position from small to large, within each chromosome)
gff <- gff[order(gff$seqname,gff$start),]

n.gene     <- nrow(gff)
gene.test2 <- rep(0,n.gene)
for (chr in unique(gff$seqname)) {
  rows.chr.i <- which(gff$seqname==chr)
  rows.chr.i <- rows.chr.i[2:(length(rows.chr.i)-1)]
  #for (j in rows.chr.i) {gene.test[j] <- (gff$end[j-1]>=gff$start[j]) & (gff$start[j+1]<=gff$end[j])}
  for (j in rows.chr.i) {gene.test2[j] <- (gff$end[j-1]>=gff$start[j+1])}
 }
gff        <- gff[which(gene.test2==0),]

### sort gff (with starting position from small to large, within each chromosome)
gff <- gff[order(gff$seqname,gff$start),]

# which(gff$start!=gff[order(gff$seqname,gff$start),][,4])
# which(gff$end<gff$start)

### define number of chromosomes, their lengths, starting positions,...
gc()
gene.lengths     <- as.integer(table(gff$seqname)) # number of genes per chromosome
gene.pos         <- c(0,cumsum(gene.lengths)[1:(nchr-1)])

# using strand information :
gff.new <- gff
for (i in 1:nchr) {
  i.index    <- gene.pos[i] + 1:(gene.lengths[i]-1)
  strand.min <- which(gff$strand[i.index]=="-")
  strand.pos <- which(gff$strand[i.index+1]=="+")
  new.end    <- apply(rbind(gff$start[i.index+1][strand.min] - 1,gff.new$end[i.index][strand.min]),2,FUN=max)
  new.end    <- apply(rbind(new.end,gff.new$end[i.index][strand.min]+max.shift),2,FUN=min)
  new.start  <- apply(rbind(gff$end[i.index][strand.pos] + 1,gff.new$start[i.index+1][strand.pos]),2,FUN=min)
  new.start  <- apply(rbind(new.start,gff.new$start[i.index+1][strand.pos]-max.shift),2,FUN=max)
  gff.new$end[i.index][strand.min]      <-  new.end
  gff.new$start[i.index+1][strand.pos]  <-  new.start
}

# which(gff.new$start!=gff.new[order(gff.new$seqname,gff.new$start),]$start)
# which(gff.new$end<gff.new$start)
### sort gff (with starting position from small to large, within each chromosome)
#gff.new <- gff.new[order(gff.new$seqname,gff.new$start),]

gene.info       <- data.frame(gene1=rep(NA,nrow(snp)),gene2=rep(NA,nrow(snp)))
#
for (i in seq(1,nrow(gff.new),by=2)) {
  gene.info$gene1[(snp$chromosome==gff.new$seqname[i]) & (snp$position %in% ((gff.new$start[i]):(gff.new$end[i])))]  <- gff.new$Name[i]
}
for (i in seq(2,nrow(gff.new),by=2)) {
  gene.info$gene2[(snp$chromosome==gff.new$seqname[i]) & (snp$position %in% ((gff.new$start[i]):(gff.new$end[i])))]  <- gff.new$Name[i]
}
# 'allign' the gene1 and gene2 columns: if there is exactly one NA for a certain gene-name
# (i.e. row in the gene.info data-frame), it should always be under gene2
na.12.ind <- (!is.na(gene.info$gene2) & is.na(gene.info$gene1))
gene.info$gene1[na.12.ind] <- gene.info$gene2[na.12.ind]
gene.info$gene2[na.12.ind] <- NA

### test
#which(apply(gene.info,1,function(x){sum(is.na(x))})==2)
#sum(apply(gene.info,1,function(x){sum(is.na(x))})==2)
################
# construction without loops ?
#gff.odd <- as.data.frame(t(gff.new[seq(1,nrow(gff.new),by=2),c(1,4:5])))
#names(gff.odd) <- gff.new$Name[seq(1,nrow(gff.new),by=2)]
#snps.odd <- lapply(gff.odd,FUN=function(x){(snp$chromosome==x[3]) & (snp$position %in% ((x[2]):(x[3])))})
##################

n.gene      <- length(unique(gff.new$Name))
genes       <- data.frame(gene.name=as.character(gff.new$Name),first.marker=rep(0,n.gene),last.marker=rep(0,n.gene),gene.length=rep(0,n.gene))
for (i in 1:nrow(genes)) {
  gene.i <- which((snp$chromosome==gff.new$seqname[i]) & (snp$position %in% ((gff.new$start[i]):(gff.new$end[i]))))
  if (length(gene.i)>0) {
    genes$first.marker[i] <- min(gene.i)
    genes$last.marker[i]  <- max(gene.i)
    genes$gene.length[i]  <- genes$last.marker[i]-genes$first.marker[i]+1
  }
}
genes$gene.name   <-  as.character(genes$gene.name)
genes <- genes[which(genes$gene.length>0),]
####test
#gene.list1 <- unique(c(gene.info$gene1[!is.na(gene.info$gene1)],gene.info$gene2[!is.na(gene.info$gene2)]))
#length(gene.list1)
#gene.list2 <- unique(genes$gene.name) #gene.list2 <- unique(gff.new$Name)
#setdiff(gene.list1,gene.list2)
#setdiff(gene.list2,gene.list1)
######

#########
n.g <- nrow(genes)
rm(n.gene,gff,gff.new)
#save(genes,n.g,gene.info,snp,n,N,plant.names,maf,file=geno.data)
save(genes,n.g,gene.info,snp,n,N,plant.names,maf,file=geno.data)

}  # end of ReadMarkerAndGeneData function


#ReadMarkerAndGeneData(gffFileName="TAIR9_GFF3_genes.gff",geno.data="wur_data.RData",snp.file.name="hapmapsnpfile_without_header.csv",header.name="hapmapsnpfile_header.csv",n.snp=214555,snp.start.col=3,ref.nr=18,maf=0)
#gffFileName="TAIR9_GFF3_genes.gff";geno.data="wur_data.RData";snp.file.name="hapmapsnpfile_without_header.csv";header.name="hapmapsnpfile_header.csv";n.snp=214555;snp.start.col=3;ref.nr=18;maf=0

# Future options / to do :
# - read kinship-matrix from a file
# - test the first lines (reading markers from a file)
#
MakeGwasObject <- function(marker.object,description="test",kinship=matrix(0),map=data.frame(),
                           pheno=data.frame(),gene.info=character(),kinship.type="ibs",
                           accession.name.type=1,allele.name.type=1,marker.name.type=2,maf=0,
                           gene.dataframe=data.frame(),generate.csv=FALSE,generate.bin=FALSE,
                           csvName=paste(description,".csv",sep=""),binName=paste(description,".bin",sep=""),
                           RimageName=paste(description,".RData",sep=""),returnObject=TRUE) {

# OBJECTIVE : given several internal and/or external objects, create a GWAS object
#
# INPUT :
#' @param marker.object :  a string, specifying a csv-file name with the marker-data
#   (row names: marker names; column-names: genotype names)
#    Alternatively, marker.object may be an R-data.frame similar to that
#  * the map (2 cols: chromosome and position) is either contained in the data-frame "marker.object",
#    or in a separate  data.frame "map" (in which case "marker.object" ONLY contains marker scores).
#    The column names should be "chromosome" and "position"; these should be the first 2 columns.
#    The positions are in base-pair or morgan. They should not be 'cumulative' over the chromosomes
#' @param  description: name of the data-set
#' @param  pheno : an R-data-frame with the phenotypic data. Importing pheno-data from a text/csv file can be done
#         using the AddPhenoData function
#' @param  accession.name.type : the type of accession names used in marker.object:
#    1= "stock.new", 2= "stock.old", 3= "ecotype.standard", 4="ecotype.yi", 5="array.id"
#' @param  allele.name.type   : the allele-encoding used in marker.object:
#    1="columbia as one",2="minor allele as one",3="other"
#' @param  marker.name.type   : the type of marker names used in marker.object:
#    1="chromosome, position and gene",2="numbers",3="other"
#' @param maf : minor allele frequency to which the markers in marker.object have been restricted
#          (in particular : before, when using the ReadMarkerAndGeneData function)
#
# * csvName
# * binName
# * RimageName : name of the R-image the GWASobject will be saved in
# * returnObject : if true, the GWASobject is returned, otherwise it is only written to the RData file
#
# OUTPUT :
# * a GWAS object (if returnObject=TRUE), where many entries and filenames will start, with "description "
# *
#
################################

# Define various factors
accession.levels   <- c("stock.new","stock.old","ecotype.standard","ecotype.yi","array.id")
allele.levels      <- c("columbia as one","minor allele as one","other")
marker.name.levels <- c("chromosome, position and gene","numbers","other")

if (is.character(marker.object)) {
  marker.object <- read.table(marker.object)
  marker.object <- as.data.frame(marker.object)
} else { # then marker.object must be an R-object. Check if it has column and row names. If not, give default names
  if (length(row.names(marker.object))==0) {row.names(marker.object)  <- paste("ind",as.character(1:nrow(marker.object)),sep="")}
  if (length(names(marker.object))==0) {names(marker.object)  <- paste("m",as.character(1:ncol(marker.object)),sep="")}
  marker.object <- as.data.frame(marker.object)
}

if (all(c("chromosome","position") %in% names(marker.object))) { # test: correct with and without map?
    map <- data.frame(chromosome=marker.object$chromosome,position=marker.object$position)
    marker.object <- marker.object[,-c(1,2)]
} else {
    if (ncol(map)>1) {names(map)[1:2]  <- c("chromosome","position")} # if the map-data-frame exists, make sure that the the first
                                                                      # 2 columns are named "chromosome" and "position"
    if (all(c("chromosome","position") %in% names(marker.object)))  {
      map <- data.frame(chromosome=map$chromosome,position=map$position)
    } else {
      map <- data.frame(chromosome=rep(1,nrow(marker.object)),position=1:nrow(marker.object))
    }
}
names(map)  <- c("chromosome","position")

###

nchr            <- length(unique(map$chromosome))
chromosomes     <- unique(map$chromosome)
chr.lengths     <- as.numeric(table(map$chromosome))
if (nchr > 1) {
  chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])
  chr.lengths.bp  <- c(0,map$position[cumsum(chr.lengths)[1:(nchr-1)]])
  cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
} else {
  cumpositions    <- map$position
  chr.pos         <- 0
  chr.lengths.bp  <- nrow(marker.object)
}
plant.names <- names(marker.object)
n           <- ncol(marker.object)
N           <- nrow(marker.object)

### principal components, as in Patterson et al 2006

gc()
PCAmatrix  <- ComputePcaMatrix(marker.object=marker.object,maf=maf)
# requiring too much memory :
# PCAmatrix  <- ComputePcaMatrix2(marker.object=marker.object,maf=maf)
PCAs       <- as.data.frame(ComputePcas(PCAmatrix)$pcas)
row.names(PCAs) <- plant.names
names(PCAs)<- paste("pca",as.character(1:ncol(PCAs)),sep="")
#PCAs <- data.frame()
### kinship matrix

if (ncol(kinship)==n) {
  A   <- kinship
} else {
    if (kinship.type=="ibs") {
        A  <- matrix(0,n,n)
        blocks <- DefineBlocks(1:N,5000)
        for (bl in 1:length(blocks)) {
          #block       <- chr.pos[CR] + 1:chr.lengths[CR]
          A  <- A + IBS(marker.object[blocks[[bl]],],normalization=1)
        }
        A  <- A/N
    }
    if (kinship.type=="id")  {
        A  <- diag(n)
    }
}
K.name          <- paste(description,"_","kinship.csv",sep="")
write.table(A, file=K.name,quote = FALSE, sep = ",", row.names = FALSE,col.names = plant.names)
AINV            <- MakeKinshipAsreml(A,genotype.names=plant.names)

### create mapframe

if (length(gene.info)>0) {
    map.frame <- cbind(map=map,data.frame(cum.position=cumpositions),gene1=gene.info$gene1,gene2=gene.info$gene2)
} else {
    map.frame <- cbind(map=map,data.frame(cum.position=cumpositions))
}

###

csv.name  <- csvName
bin.name  <- binName

if (generate.csv) {
    csv.name=csvName
    MakeCsv(markers=marker.object,plant.names=plant.names,file.name=csv.name)
} else {
    csv.name=""
}
if (generate.bin) {
    bin.name=binName
    #MakeBin(kinship.file=K.name,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
    #MakeBin(kinship=K.name,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
    MakeBin(kinship=A,plant.names=plant.names,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
} else {
    bin.name=""
}

GWAS.obj  <- list(description=description,markers=marker.object,pheno=pheno,map=map.frame,
                  kinship=A,kinship.asreml=AINV,genes=gene.dataframe,plant.names=plant.names,N=N,
                  n=n,nchr=nchr,chromosomes=chromosomes,chr.lengths.bp=chr.lengths.bp,
                  external=list(kinship.name=K.name,csv.name=csv.name,bin.name=bin.name),pca=PCAs,
                  markerInfo=list(accession.name.type=factor(x = accession.levels[accession.name.type], levels=accession.levels),
                                  allele.name.type=factor(x = allele.levels[allele.name.type], levels=allele.levels),
                                  marker.name.type=factor(x = marker.name.levels[marker.name.type], levels=marker.name.levels),maf=maf),
                  real.effects=list(locations=data.frame(NULL),sizes=data.frame(NULL))
                  )

names(GWAS.obj$map)[1:2] <- c("chromosome","position")

if (RimageName!="") {save(GWAS.obj,file=RimageName)}

if (returnObject) {return(GWAS.obj)}

}

MakeGwasObject2 <- function(geno.data="atwell_data.RData",description="test",maf=0) {
#                           csvName=paste(description,".csv",sep=""),binName=paste(description,".bin",sep=""),RimageName=paste(description,".RData",sep="")) {
# OBJECTIVE :
# load the objects created using the function ReadMarkerAndGeneData; then call
# the function MakeGwasObject to make a GWAS object
#
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *
load(geno.data)

MakeGwasObject(marker.object=snp,description=description,gene.info=gene.info,gene.dataframe=genes,
               generate.csv=T,generate.bin=T,csvName=paste(description,".csv",sep=""),
               binName=paste(description,".bin",sep=""),maf=maf,
               RimageName=paste(description,".RData",sep=""),returnObject=F)
}


##########################

#test.GWAS.obj <- make.GWAS.obj(marker.object=snp,description="test",generate.csv=T,generate.bin=T,csvName="test.csv",binName="test.bin",RimageName="test.RData")
#marker.object=snp;description="test";generate.csv=T;generate.bin=T;csvName="test.csv";binName="test.bin";RimageName="test.RData"
#data.path       <- "D:/willem/statistical.genetics.large.files/arabidopsis.data/"
#script.path     <- "D:/willem/Dropbox/research/STATISTICAL.GENETICS/arabidopsis.project/version.1.0/"
#GWAS.obj <- add.pheno.data(gwas.obj=GWAS.obj,csv.file.name="2ndHAPMAPscr.csv",add.var.means=T,mean.cols=3:4,make.pheno.image=T,pheno.image.name="charles2test.RData")
#source(paste(script.path,"functions.R",sep=""))
#setwd(data.path)

##########################
# old argument: new.gwas.image="gwas.RData"
#
# Apart from a few more comments, this new version isn't actually much different.
# Only difference: the row.names are now the accession-names + _1, _2 etc for the replicates, instead of a,b,c....

AddPhenoData  <- function(gwas.obj,csv.file.name,add.var.means=FALSE,mean.cols=0,make.pheno.image=F,pheno.image.name="pheno.RData",
                          add.normal.transform=F,which.columns.as.factor=integer(0)) {
# assumes that the working directory is correct, and that the functions script is loaded
# input  : the "gwas-object" gwas.obj and the phenotypic csv-file csv.file.name
#            which.columns.as.factor : column-numbers of the variables that are to be imported as factor
# output : the same gwas.obj, the element gwas.obj$pheno containing the data from the csv file
#           if make.pheno.image=TRUE, also an R-image containing the phenotypic data is created
#           (for consistency with earlier version this data-frame is then called data1)

#gwas.obj=GWAS.obj;add.var.means=FALSE;mean.cols=0;make.pheno.image=F;pheno.image.name="pheno.RData";add.normal.transform=F;which.columns.as.factor=integer(0)

# gwas.obj=gwas.obj;csv.file.name=csv.file.name;add.var.means=F;make.pheno.image=F;pheno.image.name=""

  if (!add.var.means) {mean.cols<-0}
  if (sum(mean.cols)==0) {add.var.means<-FALSE}
  #
  data0                   <- read.table(file=csv.file.name,sep=",",na.strings="",header=T)
  col.classes.vector      <- c("character",rep("numeric",ncol(data0)-1))
  if (length(which.columns.as.factor)>0) {col.classes.vector[which.columns.as.factor] <- "factor"}
  data0                   <- read.table(file=csv.file.name,sep=",",na.strings="",header=T,colClasses=col.classes.vector) # ,colClasses=rep("numeric",12)
  # complete the genotype-names:
  for (i in 2:nrow(data0)) {if (is.na(data0[i,1])) {data0[i,1] <- data0[i-1,1]}}
  names(data0)[1]         <- "genotype"
  # convert the accession names from factor to character
  data0$genotype          <- as.character(data0$genotype)
  data1                   <- data.frame(data0, stringsAsFactors = F)
  # remove the accessions which are not in plant.names
  data1                   <- data1[setdiff(1:nrow(data1),which(data1$genotype %in% setdiff(unique(data1$genotype),gwas.obj$plant.names))),]
  # count the number of replicates in the remaining genotypes (regardless of missing values in the phenotype)
  #rep.vec                 <- tabulate(factor(data1$genotype))
  rep.vec                 <- as.numeric(table(data1$genotype)[match(unique(data1$genotype),sort(unique(data1$genotype)))])

  # Add _1,_2, etc to the row-names # previously : add letters a,b,c... to the row-names
  ##!##
  #row.names(data1)[1:(rep.vec[1])]<- paste(data1$genotype[1:(rep.vec[1])], letters[1:(rep.vec[1])],sep="")
  row.names(data1)[1:(rep.vec[1])]<- paste(data1$genotype[1:(rep.vec[1])],"_",as.character(1:(rep.vec[1])),sep="")
  for (i in 2:length(rep.vec)) {row.names(data1)[sum(rep.vec[1:(i-1)])+ 1:(rep.vec[i])]<- paste(data1$genotype[sum(rep.vec[1:(i-1)])+ 1:(rep.vec[i])],"_",as.character(1:(rep.vec[i])),sep="")}
  # Now the row names are the accession names, where each accession name is replicated, with lower case letters as extension,
  # e.g. CS76113a,CS76113b,CS76113c

  data2                   <- data1
  # Now add extra accessions, i.e. the ones that do not occur in data1/data2, but DO occur in plant.names(2?)
  n.extra <- length(setdiff(gwas.obj$plant.names,unique(data2$genotype)))
  if (n.extra>0) {
    extra.acc               <- data.frame(matrix(NA,n.extra,ncol(data2)))
    extra.acc[,1]           <- setdiff(gwas.obj$plant.names,unique(data2$genotype))
    row.names(extra.acc)    <- paste(extra.acc[,1],"_1",sep="")
    names(extra.acc)        <- names(data2)
    data2                   <- rbind(data2,extra.acc)
    rep.vec                 <- c(rep.vec,rep(1,n.extra))
  }

  ###################################
  # Now put everything in the same order as in plant.names:

  rep.vec2     <- rep.vec[match(gwas.obj$plant.names,unique(data2$genotype))]
  plant.names2 <- rep(gwas.obj$plant.names,times=rep.vec2)
  plant.names2[1:(rep.vec2[1])] <- paste(plant.names2[1:(rep.vec2[1])], "_",as.character(1:(rep.vec2[1])),sep="")

  for (i in 2:length(rep.vec2)) {
    plant.names2[sum(rep.vec2[1:(i-1)])+ 1:(rep.vec2[i])] <-
    paste(plant.names2[sum(rep.vec2[1:(i-1)])+ 1:(rep.vec2[i])],"_",as.character(1:(rep.vec2[i])),sep="")
  }

  data3   <- data2[match(plant.names2,row.names(data2)),]

  data1   <- data3
  rep.vec <- rep.vec2

  ###################################

  if (add.var.means) {
    NC    <- ncol(data1)
    data1 <- AddMeans(input.frame=data1,col.select=mean.cols)   # the averages of the columns given in col.select are now added as extra columns
    if (add.normal.transform) {
      for (i in c(mean.cols,NC+1:length(mean.cols))) {
        data1 <- cbind(data1,qqnorm(data1[,i],plot.it=F)$x)
        names(data1)[ncol(data1)] <- paste(names(data1)[i],"_transformed",sep="")
      }
    }
  }


  if (make.pheno.image) {save(data1,rep.vec,file=pheno.image.name)}


  if (length(which.columns.as.factor)>0) {
    for (i in which.columns.as.factor) {
      data1[is.na(data1[,i]),i] <- levels(data1[,i])[1]
    }
  }

  gwas.obj$pheno  <- data1
  # old code :
  #save(gwas.obj,file=new.gwas.image)
  #assign(gwas.obj$pheno, value=data1, envir = .GlobalEnv)

return(gwas.obj)
}



#memory.size()
#gc()
#memory.size()
#gwas.obj=GWAS.obj

scan_GLS  <- function(gwas.obj,input.pheno,varcomp.file,output.file,maf=0,covariate.file="") {
# cov.string  is the string added to the scan_GLS command
#
  cov.string      <- ""
  if (covariate.file!="") {cov.string  <- paste("-fixed",covariate.file)}
  if (file.exists(gwas.obj$external$bin.name)) {
    command.string      <- paste("scan_GLS",gwas.obj$external$bin.name,input.pheno,varcomp.file,output.file,cov.string)
  } else {
    if (!file.exists(gwas.obj$external$csv.name)) {
      write.csv(gwas.obj$markers,file=gwas.obj$external$csv.name,quote=F)
    }
    gwas.obj$external$bin.name  <- paste(substr(gwas.obj$external$csv.name,1,nchar(gwas.obj$external$csv.name)-4),".bin",sep="")
    command.string      <- paste("scan_GLS",gwas.obj$external$csv.name,input.pheno,
                                 varcomp.file,output.file,cov.string,"-writebin",gwas.obj$external$bin.name)
  }

  if (maf > 0) {command.string      <- paste(command.string,"-Tminor",as.character(maf))}

  cat(command.string,"\n")
  system(command.string)
  gwa.result <- ReadGwaResult(gwas.obj=gwas.obj,output.file=output.file)

return(gwa.result)
}


# gwa.result <- ReadGwaResult(gwas.obj=GWAS.obj,output.file="D:/willem/statistical_genetics_large_files/arabidopsis_data/LFNdata/results/FT.output1_8.txt")
ReadGwaResult <- function(gwas.obj,output.file) {
  GWA.result          <- read.table(file=output.file,header=T,na.strings ="*")
  names(GWA.result) <- c("marker","Nalleles","major_allele","minor_allele","minorfreq","stat","pvalue")
  GWA.result      <- cbind(GWA.result,gwas.obj$map$chromosome,gwas.obj$map$position)
  names(GWA.result)[8:9] <- c("chromosome","position")
    if (ncol(gwas.obj$genes)!=0) {
        GWA.result          <- cbind(GWA.result,gwas.obj$map$gene1,gwas.obj$map$gene2)
        names(GWA.result)[10:11]   <- c("gene1","gene2")
    }
  GWA.result$pvalue[GWA.result$pvalue==-1]  <- 1
return(GWA.result)
}

# for scan_GLS prior to version 1.17 :
scan_GLS_sparse  <- function(gwas.obj,input.pheno,varcomp.file,output.file,covariate.file="") {
  # cov.string  is the string added to the scan.GLS command
  #
  cov.string      <- ""
  if (covariate.file!="") {cov.string  <- paste("-fixed",covariate.file)}
  if (file.exists(gwas.obj$external$bin.name)) {
    command.string      <- paste("scan_GLS_sparse",gwas.obj$external$bin.name,input.pheno,
                                 gwas.obj$external$kinship.name,varcomp.file,output.file,cov.string)
    cat(command.string,"\n")
    system(command.string, intern = TRUE, ignore.stderr = FALSE,wait = TRUE, input = NULL)
  } else {
    gwas.obj$external$bin.name  <- paste(substr(gwas.obj$external$csv.name,1,nchar(gwas.obj$external$csv.name)-4),".bin",sep="")
    command.string      <- paste("scan_GLS_sparse",gwas.obj$external$csv.name,input.pheno,gwas.obj$external$kinship.name,
                                 varcomp.file,output.file,cov.string,"-writebin",gwas.obj$external$bin.name)
    cat(command.string,"\n")
    system(command.string, intern = TRUE, ignore.stderr = FALSE,wait = TRUE, input = NULL)
  }

  gwa.result <- ReadGwaResultSparse(gwas.obj=gwas.obj,output.file=output.file)
return(gwa.result)
}

# for scan_GLS prior to version 1.17 :
ReadGwaResultSparse <- function(gwas.obj,output.file) {
    GWA.result          <- read.table(file=output.file)
    if (ncol(gwas.obj$genes)!=0) {
        GWA.result          <- cbind(GWA.result,gwas.obj$map$gene1,gwas.obj$map$gene2)
        names(GWA.result)   <- c("marker","stat","pvalue","gene1","gene2")
    } else {
        names(GWA.result)   <- c("marker","stat","pvalue")
    }
    GWA.result[GWA.result[,3]==-1,3]  <- 1
return(GWA.result)
}

DrawSubsamples <- function(pheno.dataframe,trn, plant.names,ns=100,subsample.acc=T) {
# OBJECTIVE :
# Drawing pairs of subsamples from an association panel, of half the sample-size
# (missing values not being counted; hence it matters which trait you take)
#
# INPUT :
#' @param pheno.dataframe : a phenotypic data-frame with a column named "genotype", possibly with replicates
#' @param trn             : the trait for which we need to subsample. trn is the column number in pheno.dataframe
#' @param plant.names     : the list of genotype-names
#' @param ns              : total number of subsamples. Should be even; there will be ns/2 pairs of complementary subsamples
#' @param subsample.acc   : if TRUE, subsample the genotypes (i.e. groups of replicates); otherwise
#    the individual observations are subsampled
#
# OUTPUT :
# * sub.samples           : a matrix of ns rows, one for each subsample. Every row starts with the observation numbers
#                           i.e. row numbers in pheno.dataframe) of the individuals in the subsample,
#                           then possibly followed by zeros.
# * sub.sample.size       : floor(acc.sample.size/2) (when subsample.acc=T), or
#                           floor(ind.sample.size/2) (when subsample.acc=F)
# * ind.population        : the row numbers (in pheno.dataframe) of the non-missing observations
# * ind.sample.size       : the total number of non-missing observations
# * number.of.nonmissing  : vector corresponding to plant.names.
#                           number.of.nonmissing[i] is the number of non-missing observations for genotype plant.names[i]
# * acc.population        : numbers (i.e. which entry in the vector plant.names) of genotypes
#                           with at least one non-missing observation
# * acc.sample.size       : total number of genotypes with at least one non-missing observation
############################

  ind.population  <- which(!is.na(as.numeric(pheno.dataframe[,trn]))) # the row numbers (in pheno.dataframe) of the non-missing observations
  ind.sample.size <- sum(!is.na(as.numeric(pheno.dataframe[,trn])))   # the total number of non-missing observations
  number.of.nonmissing  <- aggregate(pheno.dataframe[,trn],by=list(ordered(pheno.dataframe$genotype)),
                                     FUN=function(x){sum(!is.na(x))})[match(plant.names,sort(plant.names)),2]
  acc.population  <- which(number.of.nonmissing>0)    # numbers (i.e. which entry in the vector plant.names) of genotypes
                                                      #  with at least one non-missing observation
  acc.sample.size <- length(acc.population)           # total number of genotypes with at least one non-missing observation

  max.sub.sample.size <- max( sum(sort(number.of.nonmissing,decreasing = T)[1:floor(acc.sample.size/2)]),floor(ind.sample.size/2))
                                                      # is the maximal number of individual plants that can be contained in a subsample
  sub.samples         <- matrix(0,ns,max.sub.sample.size)

  if (subsample.acc) {
      sub.sample.size                       <- floor(acc.sample.size/2)
      for (ss in 1:(ns/2)) {
          acc.subsample1  <- sort(sample(acc.population,size=sub.sample.size))
          acc.subsample2  <- sort(sample(setdiff(acc.population,acc.subsample1),size=sub.sample.size))
          ss.size1        <- sum((pheno.dataframe$genotype %in% plant.names[acc.subsample1]) & !is.na(pheno.dataframe[,trn]))
          ss.size2        <- sum((pheno.dataframe$genotype %in% plant.names[acc.subsample2]) & !is.na(pheno.dataframe[,trn]))
          sub.samples[2*ss-1,1:ss.size1] <-  which((pheno.dataframe$genotype %in% plant.names[acc.subsample1]) & !is.na(pheno.dataframe[,trn]))
          sub.samples[2*ss,1:ss.size2]   <-  which((pheno.dataframe$genotype %in% plant.names[acc.subsample2]) & !is.na(pheno.dataframe[,trn]))
      }
  } else {
      sub.sample.size                       <-  floor(ind.sample.size/2)
      for (ss in 1:(ns/2)) {
          sub.samples[2*ss-1,1:sub.sample.size] <-  sort(sample(ind.population,size=sub.sample.size))
          sub.samples[2*ss,1:sub.sample.size]   <-  sort(sample(setdiff(ind.population,sub.samples[2*ss-1,]),size=sub.sample.size))
      }
  }
return(list(sub.samples=sub.samples,sub.sample.size=sub.sample.size,ind.population=ind.population,
       ind.sample.size=ind.sample.size,number.of.nonmissing=number.of.nonmissing,
       acc.population=acc.population,acc.sample.size=acc.sample.size))
}

# To do: rainbow colored dots for effects size
MakeLodPlot <- function(xvalues,yvalues,file.name="",jpeg.plot=T,x.lab="base pairs",
                        y.lab="-10Log(p)",x.sig=integer(0),x.effects=integer(0),
                        effects.size=numeric(0),chr.boundaries=c(0),y.thr=0) {
# OBJECTIVE : Make a plot of the LOD-profile, on screen or in a file (pdf or jpeg)
# Significant markers can be highlighted with red dots
#
# INPUT :
#' @param file.name  : if "", the plot is made to the screen. Otherwise to that file (should be .jpeg or .pdf)
#' @param jpeg.plot  : if TRUE, jpeg is produced, otherwise pdf
#' @param xvalues    : vector of cumulative marker positions
#' @param yvalues    : vector of LOD-scores
#' @param x.lab      : x-axis label
#' @param y.lab      : y-axis label
#' @param x.sig      : vector of integers, indicating which components in the vectors xvalues and yvalues are significant
#' @param x.effects  : vector of integers, indicating which components in the vector xvalues correspond to a real (known) effect
#' @param effects.size: vector of reals indicating the effect-sizes corresponding to x.effects
#' @param chr.boundaries : vector of chromosome boundaries, i.e. x-values on the same scale as xvalues
#' @param y.thr      : LOD-threshold
#
# OUTPUT : a plot
# * markers declared significant get a red dot
# * markers with a true effect get a blue dot
# * if both the real effects and the "declared significant" are given, the markers that are in the
#   intersection (i.e. true positives) get a pink dot
if (file.name!="") {
  if (jpeg.plot) {jpeg(file.name)} else {pdf(file.name)}
}
plot(x=xvalues,y=yvalues,xlab=x.lab,ylab=y.lab,type="l",lwd=0.4)
if (sum(chr.boundaries)!=0) {
  for (chr.b in chr.boundaries) {
    lines(x=rep(chr.b,2),y=c(0,max(yvalues)),col="green",lwd=1)
  }
}

class(x.sig) <- "integer"
class(x.effects) <- "integer"

if (sum(x.sig)!=0 | sum(x.effects)!=0) {
  if (sum(x.sig)!=0 & sum(x.effects)==0) {
    points(x=xvalues[x.sig],y=yvalues[x.sig],pch=20,col="red")
  }
  if (sum(x.sig)==0 & sum(x.effects)!=0) {
    points(x=xvalues[x.effects],y=yvalues[x.effects],pch=20,col="blue",lwd=2)
  }
  if (sum(x.sig)!=0 & sum(x.effects)!=0) {
    false.neg <- setdiff(x.effects,x.sig)
    false.pos <- setdiff(x.sig,x.effects)
    true.pos  <- intersect(x.sig,x.effects)
    points(x=xvalues[false.pos],y=yvalues[false.pos],pch=20,col="red")
    points(x=xvalues[false.neg],y=yvalues[false.neg],pch=20,col="blue")
    points(x=xvalues[true.pos],y=yvalues[true.pos],pch=20,col="purple",lwd=2)
  }
}
if (y.thr!=0) {lines(x=c(min(xvalues),max(xvalues)),y=rep(y.thr,2))}
if (file.name!="") {dev.off()}

} #   END OF THE FUNCTION



MakeLodPlotWithChromosomeColors <- function(xvalues,yvalues,gwas.obj,
                                           file.name="",jpeg.plot=T,
                                           x.lab="Chromosomes",
                                           y.lab=expression(-log[10](p)),
                                           plot.title="",
                                           plot.type="l",
                                           x.sig=integer(0),
                                           x.effects=integer(0),
                                           col.palette = rep(c("royalblue","maroon"),50)[1:gwas.obj$nchr],
                                           effects.size=numeric(0),
                                           chr.boundaries=c(0),
                                           y.thr=0,
                                           sign.points.thickness=0.4) {

# OBJECTIVE : Make a plot of the LOD-profile, on screen or in a file (pdf or jpeg)
# Significant markers can be highlighted with red dots
#
# INPUT :
#' @param file.name  : if "", the plot is made to the screen. Otherwise to that file (should be .jpeg or .pdf)
#' @param jpeg.plot  : if TRUE, jpeg is produced, otherwise pdf
#' @param xvalues    : vector of cumulative marker positions
#' @param yvalues    : vector of LOD-scores
#' @param x.lab      : x-axis label
#' @param y.lab      : y-axis label
#' @param plot.title : title of the plot
#' @param plot.type  : lines ("l") or dots ("d" or "p")
#' @param x.sig      : vector of integers, indicating which components in the vectors xvalues and yvalues are significant
#' @param x.effects  : vector of integers, indicating which components in the vector xvalues correspond to a real (known) effect
#' @param effects.size: vector of reals indicating the effect-sizes corresponding to x.effects
#' @param chr.boundaries : vector of chromosome boundaries, i.e. x-values on the same scale as xvalues
#' @param y.thr      : LOD-threshold
#' @param sign.points.thickness : thickness of the points that are false/true positives/negatives
#
# OUTPUT : a plot
# * markers declared significant get a red dot
# * markers with a true effect get a blue dot
# * if both the real effects and the "declared significant" are given, the markers that are in the
#   intersection (i.e. true positives) get a pink dot

#xvalues=GWAS.obj$map$cum.position;yvalues=-log10(ReplaceNaByOne(GWA.result$pvalue));plot.title=trait;gwas.obj=GWAS.obj;
#                                col.palette = rep(c("royalblue","maroon"),50)[1:GWAS.obj$nchr];
#                                x.sig=x.sig;chr.boundaries=cumsum(GWAS.obj$chr.lengths.bp)[-1];y.thr=LOD.thr.plot;
#                                x.effects=x.effects;effects.size=effects.size;plot.type="d"
#jpeg.plot=T ; x.lab="Chromosomes"; y.lab=expression(-log[10](p))
#file.name=paste("plots/",trait,suffix,".jpeg",sep="");

if (!(plot.type %in% c("l","d","p")) ) {plot.type <- "l"}
if (plot.type=="p") {plot.type <- "d"}

if (file.name!="") {
  #if (jpeg.plot) {jpeg(file.name)} else {pdf(file.name)} # ,res=300,heigth=600,width=600,quality=100
  if (jpeg.plot) {jpeg(file.name,width = 720, height = 480,quality = 100)} else {pdf(file.name)}
}

#@#
xmarks <- aggregate(x=gwas.obj$map$cum.position, by=list(factor(gwas.obj$map$chromosome)), FUN=function(x){min(x) + (max(x)-min(x))/2})[,2]
x.axis.lables <- gwas.obj$chromosomes
plot(x=xvalues,y=yvalues,xlab=x.lab,ylab=y.lab,type="n",lwd=0.4,main=plot.title,xaxt='n') # ann=FALSE
axis(1, at=xmarks, labels=x.axis.lables)

if (sum(chr.boundaries)!=0) {
  for (CHR in gwas.obj$chromosomes) {
    if (plot.type=="l") {
      lines(x=xvalues[gwas.obj$map$chromosome==CHR],y=yvalues[gwas.obj$map$chromosome==CHR],type="l",lwd=0.4,col=col.palette[which(gwas.obj$chromosomes==CHR)])
    } else {
      points(x=xvalues[gwas.obj$map$chromosome==CHR],y=yvalues[gwas.obj$map$chromosome==CHR],pch=20,lwd=0.4,col=col.palette[which(gwas.obj$chromosomes==CHR)])
    }
  }
  #for (chr.b in chr.boundaries) {
  #  lines(x=rep(chr.b,2),y=c(0,max(yvalues)),lwd=1)
  #}
} else {
    if (plot.type=="l") {
      lines(x=xvalues,y=yvalues,xlab=x.lab,ylab=y.lab,type="l",lwd=0.4,col="royalblue")
    } else {
      points(x=xvalues,y=yvalues,xlab=x.lab,ylab=y.lab,pch=20,lwd=0.4,col="royalblue")
    }
}

class(x.sig) <- "integer"
class(x.effects) <- "integer"

if (sum(x.sig)!=0 | sum(x.effects)!=0) {
  if (sum(x.sig)!=0 & sum(x.effects)==0) {
    points(x=xvalues[x.sig],y=yvalues[x.sig],pch=20,col="red",lwd=sign.points.thickness)
  }
  if (sum(x.sig)==0 & sum(x.effects)!=0) {
    points(x=xvalues[x.effects],y=yvalues[x.effects],pch=20,col="red",lwd=sign.points.thickness)
  }
  if (sum(x.sig)!=0 & sum(x.effects)!=0) {
    false.neg <- setdiff(x.effects,x.sig)
    false.pos <- setdiff(x.sig,x.effects)
    true.pos  <- intersect(x.sig,x.effects)
    points(x=xvalues[false.pos],y=yvalues[false.pos],pch=20,col="orange",lwd=sign.points.thickness)
    points(x=xvalues[false.neg],y=yvalues[false.neg],pch=20,col="red",lwd=sign.points.thickness)
    points(x=xvalues[true.pos],y=yvalues[true.pos],pch=20,col="black",lwd=sign.points.thickness)
  }
}
if (y.thr!=0) {lines(x=c(min(xvalues),max(xvalues)),y=rep(y.thr,2),lty=2,lwd=0.5)}
if (file.name!="") {dev.off()}

} #   END OF THE FUNCTION
#MakeLodPlotWithChromosomeColors(xvalues=GWAS.obj$map$cum.position,yvalues=-log10(GWA.result$pvalue),plot.title=trait,gwas.obj=GWAS.obj,col.palette = c("royalblue","maroon","royalblue","maroon","royalblue"),chr.boundaries=cumsum(GWAS.obj$chr.lengths.bp)[-1],x.effects=sim.obj$effect.locations,effects.size=sim.obj$effect.sizes,plot.type="d")



############# code found on CRAN
getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

# and here is quick parser

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

############# END of code found on CRAN


PlotMcmcResults  <- function(mcmc.results,names,effects=data.frame(),file.name="plot.pdf",no_lines=5) {
# CODE BY CHRISTIAN SCHAFER (CREST / UNIVERSITE PARIS DAUPHINE)
#
# OBJECTIVE :
# make a bar-plot of posterior inclusion probablities, based on mcmc or smc output
# *
# INPUT :
# * mcmc.results : mcmc or smc output. Nruns x p matrix, where p is the number of predictors
#                   and Nruns the number that the algorithm was run.
# * names        : vector of variable (predictor) names
# * effects      : m x 2 data-frame or matrix containing the real effects, assuming that there are m effects
#                  first column : integers indicating the location of the effects (indices corresponding to the names vector)
#                  second column: size of the effect
# * file.name    : name of the pdf-file the output is written to
#
# OUTPUT :
# a pdf-file


# boxplot data from repeated runs
boxplot = apply(mcmc.results,2,FUN=quantile,probs=c(0,0.2,0.5,0.8,1))
#boxplot = t(array(boxplot,c(length(boxplot)/5,5)))

names <- as.character(names)

if (nrow(effects) > 0) effects[,2]  <- effects[,2]/max(abs(effects[,2]))

#no_lines=1
no_bars=ceiling(length(names)/no_lines)

# create PDF-file
pdf(file=file.name, height=20, width=15,paper='a4')
#pdf(file=file.name, paper='a4')
par(mfrow=c(no_lines,1), oma=c(0, 0, 0.5, 0), mar=c(2.5, 2, 0.5, 0))

# create empty vector
empty=rep(0,length(names))

for(i in 1:no_lines) {
  start= (i-1)*no_bars + 1
  end  = min(i*no_bars, length(names))

  # create empty plot
  barplot(empty[start:end], ylim=c(0, 1), axes=FALSE, xaxs='i', xlim=c(-1, 584.2/no_lines))

  # plot effects
  if (length(effects) > 0) {
      for (i in 1:dim(effects)[1]) {
        if (start <= effects[i,1] && effects[i,1] <= end) {
          empty[effects[i,1]] = 1.05
          barplot(empty[start:end], col=rgb(1,1-abs(effects[i,2]),1-abs(effects[i,2])), axes=FALSE, add=TRUE)
          barplot(empty[start:end], col='black', axes=FALSE, angle=sign(effects[i,2])*45, density=15, add=TRUE)
          empty[effects[i,1]]=0
        }
      }
  }

  # plot results
  barplot(boxplot[,start:end], ylim=c(0, 1), names=names[start:end], las=2, cex.names=0.5, cex.axis=0.75,
          axes=TRUE, col=c('skyblue','black','white','white','black'), add=TRUE)
}
dev.off()

} # END OF FUNCTION PlotMcmcResults

mapAsDataFrame  <- function(crossObject) {
# INPUT :
# a cross-object
# OUTPUT :
# the genetic-map of the cross-object, as data-frame
# first column: chromosome; 2nd column : position
        scr          <- summary(crossObject)
        mapDataFrame <- matrix(0,sum(scr$n.mar),2)
        mapDataFrame[,1]<- rep(1:scr$n.chr,times=scr$n.mar)
        for (i in 1:scr$n.chr) {
          mapDataFrame[mapDataFrame[,1]==i,2] <- pull.map(crossObject)[[i]]
          row.names(mapDataFrame)[mapDataFrame[,1]==i] <- names(pull.map(crossObject)[[i]])
        }
        mapDataFrame <- as.data.frame(mapDataFrame)
        names(mapDataFrame) <- c("chromosome","position")
return(mapDataFrame)
}



ExtractDesignMatrixFromCross  <- function(crossObject) {
# INPUT :
# a cross-object
# OUTPUT :
# a data-frame
# first column: chromosome; 2nd column : position
  require(qtl)
  regr.frame <- data.frame(simdata=crossObject$pheno)

  for (i in 1:nchr(crossObject)) {
    n.temp  <- ncol(regr.frame)
    regr.frame <- cbind(regr.frame,as.data.frame(crossObject$geno[[i]]$data))
    names(regr.frame)[-(1:n.temp)] <- names(pull.map(crossObject)[[i]])
  }
return(regr.frame)
}

ExtractPiMassInputFilesFromCross  <- function(crossObject,geno.file="geno_piMASS.txt",pheno.file="pheno_piMASS") {
# INPUT :
# a cross-object
# OUTPUT :
# geno- and phenotype text files

  regr.frame <- ExtractDesignMatrixFromCross(crossObject)

  geno.frame <- t(2*(regr.frame[,-1]-1))
  pheno.vector<- regr.frame[,1]
  pheno.frame<- qqnorm(pheno.vector,plot.it=F)$pheno.vector
  write.table(geno.frame,file=geno.file,quote=F,row.names=T,col.names=F,sep=",")
  write.table(pheno.frame,file=pheno.file,quote=F,row.names=F,col.names=F,sep=",")

}

# to do simplify, and use the functions ExtractPiMassGenoFilesFromGwasObject and
# ExtractPiMassPhenoFilesFromGwasObject defined below
ExtractPiMassInputFilesFromGwasObject  <- function(gwas.object,trait.number,covariate.numbers=0,pca.correct=T,
                                                   geno.file="geno_piMASS.txt",pheno.file="pheno_piMASS.txt") {
# INPUT :
# a gwas-object
# trait.number : column number of the trait
# covariate.numbers :
# pca.correct :
#
# OUTPUT :
# geno- and phenotype text files
  if (pca.correct | sum(covariate.numbers)!=0) {
    regr.frame   <- data.frame(geno=gwas.object$pheno$geno,pheno=gwas.object$pheno[,trait.number],covars=gwas.object$pheno[,covariate.numbers])
    # residuals=rep(NA,nrow(gwas.object$pheno))
    if (pca.correct) {
      regr.frame   <- data.frame(regr.frame,gwas.object$pca[gwas.object$pheno$geno,])
    }
  regr.form <- as.formula(paste("pheno ~ ",paste(names(regr.frame)[3:(ncol(regr.frame))],collapse="+")))
  lin.mod   <- lm(formula=regr.form,data=regr.frame)
  #regr.frame$residuals[!is.na(regr.frame$pheno)] <- lin.mod$residuals
  gwas.object$pheno[!is.na(gwas.object$pheno[,trait.number]),trait.number] <- lin.mod$residuals
  }

  pheno.vector <- aggregate(gwas.object$pheno[,trait.number], by=list(gwas.object$pheno$geno),FUN=mean,na.rm =T)[,2]
  pheno.vector[is.nan(pheno.vector)] <- NA
  pheno.frame<- qqnorm(pheno.vector,plot.it=F)$x # OR y ?????????????????
  write.table(pheno.frame,file=pheno.file,quote=F,row.names=F,col.names=F,sep=",")

  cat("",file=geno.file)
  blocks <- DefineBlocks(1:gwas.object$N,10000)
  for (i in 1:length(blocks)) {
    maj.alleles <- round(apply(gwas.object$markers[blocks[[i]],],1,FUN=mean))
    write.table(cbind(row.names(gwas.object$markers[blocks[[i]],]),1-maj.alleles,maj.alleles,
                2*gwas.object$markers[blocks[[i]],]),file=geno.file,quote=F,row.names=F,col.names=F,sep=",",append=TRUE)
  }

}
# test:
# ExtractPiMassInputFilesFromGwasObject(gwas.object=GWAS.obj,trait.number=31)



ExtractPiMassPhenoFilesFromGwasObject  <- function(gwas.object,trait.number,covariate.numbers=0,pca.correct=T,pheno.file="pheno_piMASS.txt") {
# INPUT :
# a gwas-object
# trait.number : column number of the trait
# covariate.numbers :
# pca.correct :
#
# OUTPUT :
# phenotype text files

  if (pca.correct | sum(covariate.numbers)!=0) {
    regr.frame   <- data.frame(geno=gwas.object$pheno$geno,pheno=gwas.object$pheno[,trait.number],covars=gwas.object$pheno[,covariate.numbers])
    # residuals=rep(NA,nrow(gwas.object$pheno))
    if (pca.correct) {
      regr.frame   <- data.frame(regr.frame,gwas.object$pca[gwas.object$pheno$geno,])
    }
  regr.form <- as.formula(paste("pheno ~ ",paste(names(regr.frame)[3:(ncol(regr.frame))],collapse="+")))
  lin.mod   <- lm(formula=regr.form,data=regr.frame)
  #regr.frame$residuals[!is.na(regr.frame$pheno)] <- lin.mod$residuals
  gwas.object$pheno[!is.na(gwas.object$pheno[,trait.number]),trait.number] <- lin.mod$residuals
  }

  pheno.vector <- aggregate(gwas.object$pheno[,trait.number], by=list(gwas.object$pheno$geno),FUN=mean,na.rm =T)[,2]
  pheno.vector[is.nan(pheno.vector)] <- NA
  pheno.frame<- qqnorm(pheno.vector,plot.it=F)$x # OR y ??????????????????
  write.table(pheno.frame,file=pheno.file,quote=F,row.names=F,col.names=F,sep=",")
}
# test:
# ExtractPiMassPhenoFilesFromGwasObject(gwas.object=GWAS.obj,trait.number=31)




ExtractPiMassGenoFilesFromGwasObject  <- function(gwas.object,geno.file="geno_piMASS.txt") {
# INPUT :
# a gwas-object
#
# OUTPUT :
# geno- text file

  cat("",file=geno.file)
  #gwas.object$markers <- SwitchToMinorAlleleEncoding(marker.matrix=gwas.object$markers,two.alleles=T)
  blocks <- DefineBlocks(1:gwas.object$N,10000)
  for (i in 1:length(blocks)) {
    maj.alleles <- round(apply(gwas.object$markers[blocks[[i]],],1,FUN=mean))
    write.table(cbind(row.names(gwas.object$markers[blocks[[i]],]),1-maj.alleles,maj.alleles,
                2*gwas.object$markers[blocks[[i]],]),file=geno.file,quote=F,row.names=F,col.names=F,sep=",",append=TRUE)
  }

}
# test:
# ExtractPiMassGenoFilesFromGwasObject(gwas.object=GWAS.obj)




SwitchToMinorAlleleEncoding  <- function(marker.matrix,block.size=5000,two.alleles=F) {
# INPUT :
# a marker-matrix with zeros and ones; markers in the rows; individuals in the columns
# block.size : the operation is carried out on blocks of block.size markers
# if two.alleles=T the output-matrix is multiplied by two
# OUTPUT :
# the same marker-matrix, where, for all rows with mean larger than 0.5 the zeros and ones are interchanged
  marker.matrix   <- as.matrix(marker.matrix)
  list.of.blocks  <- DefineBlocks(1:nrow(marker.matrix),block.size=block.size)
  for (i in 1:length(list.of.blocks)){
    marker.means <- apply(marker.matrix[list.of.blocks[[i]],],1,FUN=mean)
    marker.matrix[list.of.blocks[[i]][marker.means>0.5],] <- 1-marker.matrix[list.of.blocks[[i]][marker.means>0.5],]
    if (two.alleles) {marker.matrix[list.of.blocks[[i]],] <- 2*marker.matrix[list.of.blocks[[i]],]}
  }
return(marker.matrix)
}
# qwerty <- SwitchToMinorAlleleEncoding(marker.matrix=GWAS.obj$markers,block.size=5000,two.alleles=T)

#AddPcasToGwasobject <- function(gwas.obj) {
#}


ComputePcaMatrix2  <- function(marker.object,maf=0) {
  X <- as.matrix(marker.object)
  n <- ncol(marker.object)
  minor.allele.frequencies <- apply(X,1,function(x){mean(as.numeric(x))})
  X <- X[minor.allele.frequencies >= maf & minor.allele.frequencies <= (1-maf),]
  minor.allele.frequencies <- minor.allele.frequencies[minor.allele.frequencies >= maf & minor.allele.frequencies <= (1-maf)]
  stdev <- sqrt(minor.allele.frequencies*(1-minor.allele.frequencies))
  X <- X - matrix(rep(minor.allele.frequencies,n),ncol=n)
  for (j in 1:n) {X[,j] <- X[,j]/stdev} # gwas.obj$n replaced by n
  M <- t(X) %*% X / nrow(X)
return(M)
}

ComputePcaMatrix  <- function(marker.object,maf=0) {
# INPUT :
# marker.object : a p x n matrix with markerscores that are either 0 or 1
#                 (TO DO: extend to 0,1,2)
#                 n is the number of individuals; p the number of markers
# maf : minor allele frequency : all markers with (rare) allele frequency < maf are
#        not taken into account
# OUTPUT :
# x.matrix : the matrix X= (1/p) M M' , see p. 2076 of Patterson, Price and Reich (2006)
# # removed:
# #number.of.markers : the number of markers the calculation was based on
# #(if maf=0, identical to nrow(marker.object); otherwise it may be smaller)
########
  blocks <- DefineBlocks(1:nrow(marker.object),block.size=30000)
  marker.means <- rep(0,nrow(marker.object))
  for (bl in 1:length(blocks)) {
    marker.means[blocks[[bl]]] <- apply(marker.object[blocks[[bl]],],1,mean)# / ncol(marker.object)
  }

  if (sum(marker.means<maf & marker.means>(1-maf))>0) {
    marker.object <- marker.object[marker.means>=maf & marker.means<=(1-maf),]
    marker.means <- apply(marker.object, 1, sum)/ncol(marker.object)
    }
  nr           <- nrow(marker.object)
  rescaling    <- nr*marker.means*(1-marker.means)
  nX           <- matrix(rep(0.1,(ncol(marker.object))^2),ncol=ncol(marker.object))
  ind          <- (rescaling>0)
  nX[1,1]      <- sum((marker.object[ind,1]-marker.means[ind])*(marker.object[ind,1]-marker.means[ind])/rescaling[ind])
  for (i in 2:ncol(marker.object)) {
      cat(i,"\n")
      nX[i,i] <- sum((marker.object[ind,i]-marker.means[ind])*(marker.object[ind,i]-marker.means[ind])/rescaling[ind])

      #nX[i,1:(i-1)] <- sapply(1:(i-1),function(j){   sum((marker.object[,i]-marker.means)*(marker.object[,j]-marker.means)/rescaling)}
      #nX[1:(i-1),i] <- nX[i,1:(i-1)]

      for (j in 1:(i-1)) {
        tmp <- sum((marker.object[ind,i]-marker.means[ind])*(marker.object[ind,j]-marker.means[ind])/rescaling[ind])
        nX[i,j] <- tmp; nX[j,i] <- tmp
        gc()
      }
  }
return(nX) # list(x.matrix=nX,number.of.markers=length(marker.means))
}

ComputePcas <- function(scaled.cov.matrix,sign.thr=0.05) {
# INPUT :
# scaled.cov.matrix : a n x n matrix obtained using the function ComputePcaMatrix
# sign.thr          : significance threshold used when testing for the number of significant pcas
#
# OUTPUT : a list with the components
# n.pca = the number of significant pcas
# pcas  = n x n.pca matrix, containing the scores of all individuals on these pr. components
# In a addition, a plot is made of the first two pcas (TO DO: adjust the plot; save in a file)

  if (!is.installed("RMTstat")) {install.packages("RMTstat")}
  library(RMTstat)

  m           <- nrow(scaled.cov.matrix)
  ev          <- eigen(scaled.cov.matrix)
  evalues     <- ev[[1]]
  evector1    <- ev[[2]][,1]
  evector2    <- ev[[2]][,2]
  significant <- TRUE
  k           <- 1

  while (significant) {
      meff        <- length(evalues)
      nPrime      <- (m+1)*(sum(evalues))^2 / ((m-1)*sum(evalues*evalues) - (sum(evalues))^2)
      sigma2hat   <- sum(evalues) / ((nPrime)*(m-1))
      mu.mn       <- (sqrt(nPrime-1) + sqrt(meff))^2 / nPrime
      sigma.mn    <- ((sqrt(nPrime-1) + sqrt(meff)) / nPrime) * (1/sqrt(nPrime-1) + 1/sqrt(meff))^(1/3)
      L1          <- (meff)*evalues[1] / sum(evalues)
      x.TW        <- (L1-mu.mn)/sigma.mn
      pValue      <- ptw(x.TW, beta=1, lower.tail = FALSE)
      cat(x.TW,"\t",evalues[1],"\t",pValue,"\n")
      if (pValue<sign.thr) {k <- k+1; evalues <- evalues[2:meff]} else {significant <- FALSE}
      }
  ## calculate the principal components
  pcs <- scaled.cov.matrix %*% ev[[2]][,1:k]
  ## plot :
  #couleurs <- rainbow(400)[1:m]
  #qw      <- 40
  #couleurs<- rainbow(qw+1, s = 1, v = 1, start = 0, end = max(1,qw)/(qw+1))
  #plot(x=as.vector(pcs[,1]),y=as.vector(pcs[,2]),col=couleurs)

return(list(n.pca=k,pcas=pcs,eigen.values = evalues))
}

AddPcsToGwasObj <- function(gwas.obj,maf=0.01) {
# OBJECTIVE : calculate the principal components following Price et al (2006), using the ComputePcaMatrix function, and add these to the gwas.obj
  PCAmatrix  <- ComputePcaMatrix(marker.object=gwas.obj$markers,maf=maf)
  PCAs       <- as.data.frame(ComputePcas(PCAmatrix)$pcas)
  row.names(PCAs) <- gwas.obj$plant.names
  names(PCAs)<- paste("pca",as.character(1:ncol(PCAs)),sep="")
  gwas.obj$pca <- PCAs
return(gwas.obj)
}

ConvertACGTfactorLevelsTo1234 <- function(x) {
# INPUT : a vector of A,C,G and T's , of class factor
# OUTPUT : the same vector with the levels A,C,G and T replaced by 1,2,3,4
  levels(x)[levels(x)=="A"]<-"1"
  levels(x)[levels(x)=="C"]<-"2"
  levels(x)[levels(x)=="G"]<-"3"
  levels(x)[levels(x)=="T"]<-"4"

return(as.integer(as.character(x)))
}

# to do: n.sim >1
SimulatePhenoGivenGeno1 <- function(marker.object,effect.variance=1,n.effect=30,PVE=0.5,replicates=1,
                                    equal.effect.sizes=F,equal.var.explained=F,min.distance=1,n.sim=1,
                                    group.vector=numeric(0),min.effect.distance=0) {

#browser()
# OBJECTIVE : as in section 5.1 of Guan and Stephens (2011) :
# simulation of a phenotype corresponding to real genotypic GWAS data. Gaussian effects
# are simulated, and the residual variance is chosen in order to have a prespecified value
# of the PVE (explained below). Given these effects and residual variance, a normally
# distributed phenotype is simulated.
#
# INPUT :
# * marker.object :  p x n.ind matrix or data.frame of 0/1 marker-scores
# * effect.variance : variance of the (normally distributed) effects
# * n.effect : number of effects
# * PVE : proportian of explained variance; this mimics the heritability
#   see Guan and Stephens (2011) for more details
# * replicates : number of replicates for each genotype
# * n.sim : number of simulated data sets
# * equal.effect.sizes : if TRUE, effect.sizes     <- rep(sqrt(effect.variance),n.effect)
#
# OUTPUT :
# *
# *
# *
  n.ind            <- ncol(marker.object)*replicates
  if (length(group.vector)==0) {
    effect.locations <- sort(sample(1:nrow(marker.object),size=n.effect))
  } else {
    group.labels <- unique(group.vector)
    if (n.effect > length(group.labels)) {cat("ERROR: more effects than groups.")}
    selected.groups <- sort(sample(group.labels,size=n.effect))
    effect.locations <- numeric(0)
    for (j in selected.groups) {
      effect.locations <- c(effect.locations,sort(sample((1:nrow(marker.object))[group.vector==j],size=1)))
    }
  }
  if (min.effect.distance>0 & n.effect>1) {
    while ( sum((effect.locations-c(-min.effect.distance,effect.locations[1:(n.effect-1)])) < min.effect.distance)>0 ) {
      #n.effect <- sum((effect.locations-c(-min.effect.distance,effect.locations[1:(n.effect-1)])) >= min.effect.distance)
      effect.locations <-  effect.locations[(effect.locations-c(-min.effect.distance,effect.locations[1:(n.effect-1)])) >= min.effect.distance]
      n.effect <- length(effect.locations)
    }
  }
  if (n.effect>0) {
    if (equal.effect.sizes) {
      effect.sizes     <- rep(sqrt(effect.variance),n.effect)
      if (equal.var.explained) {
        effect.fr <- apply(marker.object[effect.locations,],1,mean)
        effect.sizes <- effect.sizes / sqrt(effect.fr*(1-effect.fr))
      }
    } else {
      effect.sizes     <- rnorm(sd=sqrt(effect.variance),n=n.effect)
    }
    predicted        <- rep(as.numeric(as.matrix(t(as.matrix(marker.object[effect.locations,]))) %*% matrix(effect.sizes,ncol=1)),each=replicates)
  } else {
    predicted        <- rep(0,n.ind)
    effect.sizes     <- numeric()
  }
  emp.effect.var   <- var(predicted) # (predicted)^2 / n.ind
  #res.variance     <- emp.effect.var * (1-PVE) / PVE
  if (emp.effect.var==0) {
    res.variance     <- 1
   } else {
     res.variance     <- emp.effect.var * (1-PVE)
   }
  y                <- data.frame(matrix(0,n.ind,n.sim))
  names(y)         <- paste("sim",as.character(1:n.sim),sep="")
  for (i in 1:n.sim) {y[,i] <- predicted + rnorm(sd=sqrt(res.variance),n=n.ind)}
  #y                <- predicted + rnorm(sd=sqrt(res.variance),n=n.ind)
return(list(phenotype=y,effect.locations=effect.locations,effect.sizes=effect.sizes))
}

SimulatePhenoGivenGeno2 <- function(gwas.obj,h2=0.5,n.rep=1,qtl.loc=1000,n.sim=1,pr.gen.var=0.1,subsample.size.sim=gwas.obj$n,effect.size=1) {

# simulation of a phenotype corresponding to real genotypic GWAS data.
# Following MacKenzie and Hackett (2011, appendix; assuming no additional variation due to blocks)
#
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
# * n.rep : number of replicates for each genotype
# * n.sim : number of simulated data sets
# * qtl.loc : the row in gwas.obj$markers where the qtl is
# * h2 : heritability
# * pr.gen.var : proportion of the genetic (not total) variance explained by the qtl
#
# OUTPUT :
# *
# *
# *
  #a          <- 1   # effect size of a single minor allele (so multiply this by two)  ===> WRONG : a is indeed the effect of being homozygous
  if (subsample.size.sim > gwas.obj$n) {subsample.size.sim <- gwas.obj$n}
  pheno.frame        <- data.frame(matrix(0,n.rep*gwas.obj$n,n.sim+1))
  pheno.frame[,1]    <- rep(gwas.obj$plant.names,each=n.rep)
  names(pheno.frame) <- c("genotype",paste("sim",as.character(1:n.sim),sep=""))

  fr      <- mean(as.numeric(gwas.obj$markers[qtl.loc,]))
  sigmaC2  <- fr * (1-fr) * gwas.obj$n * effect.size^2
  sigmaA2    <- ((1-pr.gen.var)/pr.gen.var) * sigmaC2 / KinshipTransformUnscaled(gwas.obj$kinship)    # polygenic variance
  sigmaG2 <- sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship)
  sigmaE2 <- (sigmaG2 + sigmaC2) * (1-h2) / (h2 * (gwas.obj$n-1)) #

  #var.tot <- sigmaG2 + sigmaE2*(gwas.obj$n-1) + sigmaC2 # a^2*fr*(1-fr)
  #var.tot
  #var.tot/(gwas.obj$n-1)
  #(sigmaG2 + sigmaC2) / var.tot
  #sigmaC2 #a^2*fr*(1-fr) * gwas.obj$n
  #sigmaG2
  #sigmaE2*(gwas.obj$n-1)

  for (trait in 1:n.sim) {
    G <- MatrixRoot(gwas.obj$kinship) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))
    z <- c(rep(c(rep(1,n.rep),rep(0,n.rep*gwas.obj$n)),gwas.obj$n-1),rep(1,n.rep))
    Z <- matrix(z,ncol=gwas.obj$n)
    G <- as.numeric(Z %*% G)
    pheno.frame[,1+trait] <- G + rnorm(n.rep * gwas.obj$n, sd=sqrt(sigmaE2)) +  effect.size*rep(as.numeric(gwas.obj$markers[qtl.loc,]),each=n.rep)
    if (subsample.size.sim!=gwas.obj$n) {pheno.frame[sample(1:nrow(pheno.frame),size=nrow(pheno.frame)-subsample.size.sim),1+trait] <- NA}

  }

return(list(pheno.frame=pheno.frame,sigmaA2=sigmaA2,sigmaG2=sigmaG2/(gwas.obj$n-1),sigmaE2=sigmaE2,effectVar=sigmaC2/gwas.obj$n))
}


SimulatePhenoGivenGeno2Old <- function(marker.object,sigma2a=1,sigma2e=1,h=0.5,mu=10,A.matrix,replicates=1,n.sim=1) {
# OBJECTIVE :
# simulation of a phenotype corresponding to real genotypic GWAS data.
# Following MacKenzie and Hackett (2011, appendix; assuming one block)
#
# INPUT :
# * marker.object : p x n.ind  matrix or data.frame of 0/1 marker-scores
# * replicates : number of replicates for each genotype
# * n.sim : number of simulated data sets
  n.ind            <- ncol(marker.object)*replicates
  n.mar            <- nrow(marker.object)
  #
  y                <- matrix(0,n.ind,n.sim)
  effect.loc       <- rep(0,n.sim)
  marker.freq      <- rep(0,n.sim)
  a2               <- rep(0,n.sim)
  #
  for (i in 1:n.sim) {
    effect.loc[i]  <- sample(1:n.mar,size=1)
    marker.freq[i] <- mean(as.numeric(marker.object[effect.loc[i],]),na.rm=T)
    a2[i]          <- (1/(n.ind*marker.freq[i]*(1-marker.freq[i]))) * ((n.ind-1)*sigma2e*h^2/(1-h^2) + sigma2a*(n.ind-(1/n.ind)*sum(A.matrix)))
    #
    predicted      <- rep(as.numeric(as.matrix(t(marker.object[effect.loc[i],])) %*% matrix(sqrt(a2[i]),ncol=1)),each=replicates)
    predicted      <- predicted + mu
    y[,i]          <- predicted + rnorm(sd=sqrt(sigma2e),n=n.ind) + rep(mvrnorm(n = 1, mu=rep(0,n.ind/replicates), Sigma=2*sigma2a*A.matrix, tol = 1e-6),each=replicates)
  # as.numeric(... %*% ... )
  }
return(list(phenotype=y,effect.locations=effect.loc,qtl.freq=marker.freq,effect.sizes=sqrt(a2)))
}

SimulatePhenoGivenGeno3 <- function(gwas.obj,h2=0.5,n.rep=1,qtl.loc=1000,n.sim=1,pr.gen.var=0.1,subsample.size.sim=gwas.obj$n) {

# simulation of a phenotype corresponding to real genotypic case-control data.
# for a population without structure
# N.B. no replicates !
#
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
# * n.sim : number of simulated data sets
# * qtl.loc : the row in gwas.obj$markers where the qtl is
# * h2 : heritability
# * pr.gen.var : proportion of the total variance explained by the qtl
#
# OUTPUT :
# *
# *
# *
  a          <- 1   # effect size of a single minor allele (so multiply this by two)
  if (subsample.size.sim > gwas.obj$n) {subsample.size.sim <- gwas.obj$n}
  pheno.frame        <- data.frame(matrix(0,gwas.obj$n,n.sim+1))
  pheno.frame[,1]    <- rep(gwas.obj$plant.names)
  names(pheno.frame) <- c("genotype",paste("sim",as.character(1:n.sim),sep=""))

  #fr      <- mean(as.numeric(gwas.obj$markers[qtl.loc,]))
  #sigmaC2  <- fr * (1-fr) * gwas.obj$n * a^2
  #sigmaE2 <- (sigmaC2) * (1-h2) / (h2 * (gwas.obj$n-1)) #
  sigmaE2  <- var(a*(as.numeric(gwas.obj$markers[qtl.loc,]))) * (1/h2 - 1)

  for (trait in 1:n.sim) {
    pheno.frame[,1+trait] <- rnorm(gwas.obj$n, sd=sqrt(sigmaE2)) +  a*(as.numeric(gwas.obj$markers[qtl.loc,]))

    threshold <- sort(pheno.frame[,1+trait],decreasing=T)[mean(gwas.obj$pheno$disease)*gwas.obj$n]
    pheno.frame[,1+trait] <- as.numeric(pheno.frame[,1+trait] > threshold)

    if (subsample.size.sim!=gwas.obj$n) {pheno.frame[sample(1:nrow(pheno.frame),size=gwas.obj$n-subsample.size.sim),1+trait] <- NA}
  }

return(list(pheno.frame=pheno.frame,sigmaE2=sigmaE2))
}


SimulatePhenoGivenGeno4 <- function(gwas.obj,n.sim=100,subsample.size.sim=round(gwas.obj$n/2)) {

# simulation of a phenotype corresponding to real genotypic case-control data for a population without structure,
# using subsampling.
# N.B. no replicates !
#
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores.
#              Its pheno data-frame should have a column 'disease', without missing values
# * n.sim : number of simulated data sets
# * subsample.size.sim : sample size
# OUTPUT :
# *
# *
# *
  #a          <- 1   # effect size of a single minor allele (so multiply this by two)
  if (subsample.size.sim > gwas.obj$n) {subsample.size.sim <- gwas.obj$n}
  pheno.frame        <- data.frame(matrix(0,gwas.obj$n,n.sim+1))
  pheno.frame[,1]    <- gwas.obj$plant.names
  names(pheno.frame) <- c("genotype",paste("sim",as.character(1:n.sim),sep=""))

  for (trait in 1:n.sim) {
    pheno.frame[,1+trait] <- gwas.obj$pheno$disease
    #if (subsample.size.sim!=gwas.obj$n) {pheno.frame[sample(1:nrow(pheno.frame),size=gwas.obj$n-subsample.size.sim),1+trait] <- NA}
    if (subsample.size.sim!=gwas.obj$n) {
      #pheno.frame[sample(1:nrow(pheno.frame),size=gwas.obj$n-subsample.size.sim),1+trait] <- NA
      n1 <- sum(pheno.frame[,1+trait])
      n0 <- nrow(pheno.frame) - n1
      n1.sample <- round(subsample.size.sim * n1 / (n0+n1))
      n0.sample <- subsample.size.sim - n1.sample
      pheno.frame[sample(which(pheno.frame[,1+trait]==1),size=n1.sample),1+trait] <- 200
      pheno.frame[sample(which(pheno.frame[,1+trait]==0),size=n0.sample),1+trait] <- 100
      pheno.frame[!(pheno.frame[,1+trait] %in% c(100,200)),1+trait] <- NA
      pheno.frame[which(pheno.frame[,1+trait]==100),1+trait] <- 0
      pheno.frame[which(pheno.frame[,1+trait]==200),1+trait] <- 1
      #pheno.frame[sample(1:nrow(pheno.frame),size=gwas.obj$n-subsample.size.sim),1+trait] <- 100
    }
  }

return(list(pheno.frame=pheno.frame))
}



DefineGenomicRegions <- function(chr.vector,block.size=100) {
  n.chr  <- length(unique(chr.vector))
  chr.nrs<- sort(unique(chr.vector))
  output.list1 <- list()
  for (i in chr.nrs) {
    j <- which(chr.nrs==i)
    n.marker <- sum(chr.vector==i)
    k1 <- ceiling(n.marker/block.size)
    difference <- round(block.size/2)
    output.list1[[j]] <- data.frame(begin=1+(0:(k1-1))*block.size,end=(1:k1)*block.size)
    output.list1[[j]][nrow(output.list1[[j]]),2] <- n.marker
    output.list1[[j]] <- output.list1[[j]] + sum(chr.vector<i)
  }
output.list1
}
# test:
#DefineGenomicRegions(chr.vector=GWAS.obj$map$chromosome,block.size=10000)

DefineGenomicRegionsWithOverlap <- function(chr.vector,block.size=100) {
  if (ceiling(block.size/2)!=block.size/2) {block.size <- block.size + 1}
  n.chr  <- length(unique(chr.vector))
  chr.nrs<- sort(unique(chr.vector))
  output.list1 <- data.frame(begin=integer(0),end=integer(0))
  for (i in chr.nrs) {
    j <- which(chr.nrs==i)
    n.marker <- sum(chr.vector==i)
    k <- ceiling(2*n.marker/block.size)
    difference <- block.size/2
    new.output.list <- data.frame(begin=1+(0:(k-1))*(block.size/2),end=block.size+(0:(k-1))*(block.size/2))
    nr <- nrow(new.output.list)
    new.output.list[(nr-1):nr,2] <- n.marker
    new.output.list <- new.output.list + sum(chr.vector<i)
    output.list1 <- rbind(output.list1,new.output.list)
  }
return(output.list1)
}

haldane <- function(x){exp(-2*x/100)}

computeVariance <- function(effectsFrame) {
# useful for heritability calculations
  x <- matrix(effectsFrame$location,ncol=nrow(effectsFrame))
  a <- matrix(outer(t(x),x,FUN="-"),ncol=nrow(effectsFrame))
  y <- matrix(effectsFrame$size,ncol=1)
  return(t(y) %*% haldane(abs(a)) %*% y)
}

MoskvinaCorrection <- function(marker.frame,number.of.replicates=rep(1,ncol(marker.frame)),
                               inv.cor.matrix=diag(rep(1,sum(number.of.replicates))),b.size=100,alpha=0.05) {
# OBJECTIVE :
# *
# *
# INPUT :
# * marker.frame : markers in the rows; genotypes in the columns. Data from one chromosome !
# * number.of.replicates : a vector of length ncol(marker.frame).
#   For every genotype, it gives the number of observations (individuals of that genotype)
# * inv.cor.matrix : the inverse of the correlation matrix of the individual observations.
#   Should have dimension sum(number.of.replicates) x sum(number.of.replicates)
# * b.size : the block-size : for every marker, we look at the correlation with the b.size preceding markers
# * alpha : desired FWE-control
# OUTPUT :
# *
# *
# *
# USES the functions MatrixRoot and DefineBlocks
N     <- nrow(marker.frame)
n.geno<- ncol(marker.frame)
Kappa <- rep(1,N)

v.blocks <- DefineBlocks(2:N,block.size=b.size)
n.blocks <- length(v.blocks)
h.blocks <- v.blocks
h.blocks <- lapply(h.blocks,FUN=function(x){ymin <- min(x); ymax<- max(x); (ymin-b.size):(ymax-1)})
h.blocks[[1]] <- 1:(b.size+1)
#h.blocks[[n.blocks]] <- (min(h.blocks[[n.blocks]])):(max(max(h.blocks[[n.blocks]]),N))

INVROOT <- MatrixRoot(inv.cor.matrix)

for (b in 1:n.blocks) {
  X <- as.matrix(marker.frame[(min(h.blocks[[b]])):(max(v.blocks[[b]])),rep(1:n.geno,times=number.of.replicates)])
  #if (doubleGeno) {X <- 2*X}
  # premultiply X with the square root of inv.cor.matrix, and transpose
  X <- INVROOT %*% t(X)
  #A <- cor(X[,v.blocks[[b]]],X[,h.blocks[[b]]])
  A <- cor(X[,(ncol(X)-length(v.blocks[[b]])+1):(ncol(X))],X[,1:length(h.blocks[[b]])]) # b.size + 1:b.size],X[,1:(b.size+1)
  if (b!=1) {
    A[1:min(b.size,length(v.blocks[[b]])),1:min(b.size,length(v.blocks[[b]]))][lower.tri(A[1:min(b.size,length(v.blocks[[b]])),1:min(b.size,length(v.blocks[[b]]))])] <- 0
    A[1:min(b.size-1,length(v.blocks[[b]])),max(1,ncol(A)-b.size+2):ncol(A)][upper.tri(A[1:min(b.size-1,length(v.blocks[[b]])),max(1,ncol(A)-b.size+2):ncol(A)],diag=T)] <- 0
  } else {
    A[!lower.tri(A,diag=T)]<-0
  }
  # detail : diag=F for last block, in upper.tri !
  A <- abs(A)
  Kappa[v.blocks[[b]]] <- apply(A,1,max)
}

Kappa <- sqrt(1-Kappa^{-1.31*log(alpha,base=10)})
return(sum(Kappa))
}

GaoCorrection <- function(marker.frame,number.of.replicates=rep(1,ncol(marker.frame)),cut.off=0.995,
                               inv.cor.matrix=diag(rep(1,sum(number.of.replicates))),doubleGeno=TRUE) {
# INPUT :
# * marker.frame : markers in the rows; genotypes in the columns. Data from one chromosome !
# * number.of.replicates : a vector of length ncol(marker.frame).
#   For every genotype, it gives the number of observations (individuals of that genotype)
# * inv.cor.matrix : the inverse of the correlation matrix of the individual observations.
#   Should have dimension sum(number.of.replicates) x sum(number.of.replicates)
# * alpha : desired FWE-control
# OUTPUT :
# *
# *
# *
# USES the functions MatrixRoot and DefineBlocks
N         <- nrow(marker.frame)
n.geno    <- ncol(marker.frame)
INVROOT   <- MatrixRoot(inv.cor.matrix)
blocksize <- 2000
matrix.blocks <- DefineBlocks(indices=1:N,block.size=blocksize)
n.blocks      <- length(matrix.blocks)
new.matrix    <- matrix(0,sum(number.of.replicates),N)

for (b in 1:n.blocks) {
  X <- t(as.matrix(marker.frame[matrix.blocks[[b]],rep(1:n.geno,times=number.of.replicates)]))
  if (doubleGeno) {X <- 2*X}
  # premultiply X with the square root of inv.cor.matrix, and transpose
  X <- INVROOT %*% X
  new.matrix[,matrix.blocks[[b]]] <- as.matrix(as.data.frame(lapply(data.frame(X),function(x){(x-mean(x))/sd(x)})))
}
sin.vals   <- svd(new.matrix)$d
eigen.vals <-  sin.vals^2
which(cumsum(eigen.vals)/sum(eigen.vals) > 0.999)
Meff <- min(which(cumsum(eigen.vals)/sum(eigen.vals) > cut.off))

return(Meff)
}

CalculateEffect <- function(gwas.obj,snp.name,trait.number) {
  effect.size <- NA
  if (!is.na(mean(gwas.obj$pheno[gwas.obj$markers[snp.name,gwas.obj$pheno$genotype]==1,trait.number],na.rm=T)) & !is.na(mean(gwas.obj$pheno[gwas.obj$markers[snp.name,gwas.obj$pheno$genotype]==0,trait.number],na.rm=T))) {
    effect.size <- mean(gwas.obj$pheno[gwas.obj$markers[snp.name,gwas.obj$pheno$genotype]==1,trait.number],na.rm=T) - mean(gwas.obj$pheno[gwas.obj$markers[snp.name,gwas.obj$pheno$genotype]==0,trait.number],na.rm=T)
  }
effect.size/2
}

ReadGenstatSimulatedCross <- function(geno.file,map.file,pheno.file,rqtl.file,trait.name="test") {
# OBJECTIVE : read  a geno, map , and pheno file (typically the DH-simulation by Marcos, done in Genstat),
#             write it to a csvr file (see www.rqtl.org) import it in R-qtl
# INPUT :
# *
# *
# *
# OUTPUT :
# * a cross-object
  require(qtl)
  map             <- read.table(map.file,sep="\t")
  map[,1]         <- as.character(map[,1])
  geno            <- read.table(geno.file,header=T,sep="\t")
  con             <- file(description=pheno.file, open = "r")
  pheno           <- read.table(con,sep=",", nrows = 1) # ,colClasses=c(rep("character",n+2))
  pheno           <- read.table(con,sep=",", nrows = dim(geno)[1] + 5)
  close(con)
  #csvr format
  cat(trait.name,"","",pheno[,2],file=rqtl.file,append=F,sep=",")
  cat("\n",file=rqtl.file,append=T,sep="")
  write.table(cbind(map,t(geno[-(1:2),-1])),quote=F,na="",sep=",",append=T,file=rqtl.file,row.names=F,col.names=F)
  sim.a <- read.cross("csvr", dir=getwd(), rqtl.file,genotypes=c("1","2"))
return(sim.a)
}

GetResiduals <- function(y,X,lm.family=1) {
# INPUT : a vector y and a design matrix or data.frame X, without intercept
#   family = 1 standard regression. family!=1 : logit
# OUTPUT : residuals of regression of y on to X
  Data       <- data.frame(cbind(y,X))#if (!is.data.frame)
  lm.formula <- paste(names(Data)[1],"~",paste(names(Data)[-1],collapse="+"))
  if (lm.family==1) {
    lm.fit     <- glm(formula=lm.formula,data=Data,family=gaussian(link = "identity"))
  } else {
    lm.fit     <- glm(formula=lm.formula,data=Data,family=binomial(link = "logit"))
  }
return(as.numeric(lm.fit$residuals))
}

NormalQuantileTransform <- function(y) {
# the normal quantile transform as in Guan and Stephens (2010)
return(qqnorm(y,plot.it=F)$x)
}

KinshipTransform <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}

KinshipTransformUnscaled <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn))
}

CompareResults <- function(effect.locations.ind,estimated.locations.ind,tolerance=0) {
# N.B. Assumes that all effects are marker positions
# N.B. The tolerance argument ignores chromosome boundaries: for each true effect, the interval of width 2*tolerance around the effect is assumed to be on the same chromosome
# N.B. assumes that 2*tolerance is smaller than the smallest distance occurring between 2 effects (in numbers of markers)
#
# OBJECTIVE : compare the locations of true effects given in effect.locations.ind with the estimated locations in estimated.locations.ind,
#             and calculate the number of false positives and negatives
# INPUT :
# * effect.locations.ind  : binary vector, its length being the number of markers. Each element is TRUE if there is a (true) effect at the corresponding marker; F otherwise
# * estimated.locations.ind : binary vector, its length being the number of markers. Each element is TRUE if, according to the estimate, there is an effect
# * tolerance :
# OUTPUT :
# the numbers of respectively false positives and false negatives
#
  nm                  <- length(effect.locations.ind)
  n.effect            <- sum(effect.locations.ind)
  effect.locations    <- which(effect.locations.ind)
  estimated.locations <- which(estimated.locations.ind)
  #
  effect.regions      <- integer(0)
  for (discr in -tolerance:tolerance) {effect.regions <- c(effect.regions,effect.locations+discr)}
  effect.regions      <- sort(effect.regions)
  effect.regions.ind  <- rep(F,nm)
  effect.regions.ind[effect.regions] <- T
  #
  region.number       <- rep(0,nm)
  region.number[effect.regions.ind]  <- rep(1:n.effect,each=2*tolerance+1)
  region.frame        <- data.frame(region.number=region.number,estimated=as.vector(estimated.locations.ind))
  #
  comparison.frame    <- aggregate(region.frame$estimated,by=list(region.number),FUN=sum)
  fp                  <- sum(comparison.frame[1,2]) + sum(pmax(rep(0,n.effect),comparison.frame[2:(1+n.effect),2]-rep(1,n.effect))) #rep(1,n.effect)
  fn                  <- 1-comparison.frame[2:(1+n.effect),2]
  fn[fn < 0]          <- 0
#
return(list(number.of.false.positives=fp,number.of.false.negatives=sum(fn)))
}

ExpandIntegers <- function(sequ) {
  result <- integer(0)
  for (i in 1:length(sequ)) {result <- c(result,1:sequ[i])}
return(result)
}

TransformGwasObject <- function(gwas.obj,marker.selection=1:gwas.obj$N,modify.external=T) { # ,ind.selection=1:gwas.obj$n
#browser()
# OBJECTIVE :
# *
# *
# INPUT :
# * gwas.obj, GWAS.object with 'markers' field being a p x n matrix or data.frame (n=number of genotypes)
# *  modify.external : if TRUE, re-create csv- and bin files of the genotypic data
# *
# OUTPUT :
# *
# *
# *
  if (nrow(gwas.obj$markers)>0) {  # arabidopsis data
    new.gwas.obj <- list()
    new.gwas.obj$description     <- gwas.obj$description
    new.gwas.obj$markers         <- gwas.obj$markers[marker.selection,] # ind.selection
    new.gwas.obj$pheno           <- gwas.obj$pheno
    new.gwas.obj$map             <- gwas.obj$map[marker.selection,]
    #new.gwas.obj$kinship         <- gwas.obj$kinship
    new.gwas.obj$kinship         <- IBS(gwas.obj$markers)
    new.gwas.obj$kinship.asreml  <- gwas.obj$kinship.asreml
    new.gwas.obj$genes           <- gwas.obj$genes
    new.gwas.obj$plant.names     <- gwas.obj$plant.names
    new.gwas.obj$N               <- length(marker.selection)
    new.gwas.obj$n               <- gwas.obj$n # length(ind.selection)
    new.gwas.obj$nchr            <- length(unique(new.gwas.obj$map$chromosome))
    new.gwas.obj$chromosomes     <- sort(unique(new.gwas.obj$map$chromosome))
    #
    new.chr.lengths     <- as.numeric(table(new.gwas.obj$map$chromosome))
    if (new.gwas.obj$nchr > 1) {
      chr.pos         <- c(0,cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)])
      chr.lengths.bp  <- c(0,new.gwas.obj$map$position[cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)]])
      #cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
    } else {
      #cumpositions    <- map$position
      #chr.pos         <- 0
      chr.lengths.bp  <- 0#nrow(marker.object)
    }
    new.gwas.obj$chr.lengths.bp     <- chr.lengths.bp
    #
    if (modify.external) {
      new.gwas.obj$external        <- gwas.obj$external
      new.gwas.obj$external$csv.name <- paste("modified_",gwas.obj$external$csv.name,sep="")
      new.gwas.obj$external$bin.name <- paste("modified_",gwas.obj$external$bin.name,sep="")
      #
      MakeCsv(markers=new.gwas.obj$markers,plant.names=new.gwas.obj$plant.names,file.name=new.gwas.obj$external$csv.name)
      MakeBin(kinship=new.gwas.obj$kinship,plant.names=gwas.obj$plant.names,pheno=gwas.obj$pheno,
              csv.file.name=new.gwas.obj$external$csv.name,bin.file.name=new.gwas.obj$external$bin.name)    # old argument : markers=new.gwas.obj$markers,
      #MakeBin(markers=new.gwas.obj$markers,kinship.file=new.gwas.obj$external$kinship.name,csv.file.name=new.gwas.obj$external$csv.name,bin.file.name=new.gwas.obj$external$bin.name)
      #
    }
    new.gwas.obj$pca             <- gwas.obj$pca

    new.gwas.obj$markerInfo      <- gwas.obj$markerInfo
    new.gwas.obj$real.effects    <- gwas.obj$real.effects
    #new.gwas.obj$     <- gwas.obj$

  } else { # WTCCC-data

    new.gwas.obj <- list()
    #new.gwas.obj$description     <- gwas.obj$description
    new.gwas.obj$pheno           <- gwas.obj$pheno
    new.gwas.obj$map             <- gwas.obj$map[marker.selection,]
    new.gwas.obj$plant.names     <- gwas.obj$plant.names
    new.gwas.obj$N               <- length(marker.selection)
    new.gwas.obj$n               <- gwas.obj$n # length(ind.selection)
    new.gwas.obj$nchr            <- length(unique(new.gwas.obj$map$chromosome))
    new.gwas.obj$chromosomes     <- sort(unique(new.gwas.obj$map$chromosome))
    #
    new.chr.lengths     <- as.numeric(table(new.gwas.obj$map$chromosome))
    if (new.gwas.obj$nchr > 1) {
      chr.pos         <- c(0,cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)])
      chr.lengths.bp  <- c(0,new.gwas.obj$map$position[cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)]])
      #cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
    } else {
      #cumpositions    <- map$position
      #chr.pos         <- 0
      chr.lengths.bp  <- 0#nrow(marker.object)
    }
    new.gwas.obj$chr.lengths.bp     <- chr.lengths.bp
    #
    new.gwas.obj$external        <- gwas.obj$external
    new.plink.name                  <- paste("modified_",gwas.obj$external$plink.name,sep="")
    new.gwas.obj$external$plink.name <- new.plink.name
    #
    #if (!all(file.exists(paste(new.plink.name,".bim",sep="")),file.exists(paste(new.plink.name,".bed",sep="")),file.exists(paste(new.plink.name,".fam",sep="")))) {
      temp.snp.file <- "tempsnpfile.txt"
      write.table(gwas.obj$map$snp.name[marker.selection],file=temp.snp.file,quote=F,col.names=F,row.names=F)
      system(paste("plink --bfile",gwas.obj$external$plink.name,"--make-bed --extract",temp.snp.file,"--out",new.plink.name))
    #}
    system(paste("plink --bfile",new.plink.name,"--recodeA --out",new.plink.name))
    new.gwas.obj$markers <- t(read.table(file=paste(new.plink.name,".raw",sep=""),header=T)[,-(1:6)])#t(read.table(file="modified_CD_and_wtc1.raw",header=T)[,-(1:6)])
    new.gwas.obj$pca             <- gwas.obj$pca
  }
  if (length(gwas.obj$groups)>0) {
    group.selection <- apply(gwas.obj$groups$position.matrix,1,FUN=function(x){sum(marker.selection %in% (x[1]):(x[2]))})
    G.temp <- gwas.obj$groups$position.matrix[group.selection>0,]
    G.temp[1,1] <- marker.selection[1]
    G.temp[nrow(G.temp),2] <- marker.selection[length(marker.selection)]
    G.temp <- G.temp - G.temp[1,1] + 1
    new.gwas.obj <- AddGroupsToGwasObj(G.temp,new.gwas.obj)
  }
return(new.gwas.obj)
}

pos.exp <- function(x,exponent=1) {x.out <- x; x.out[x.out!=0] <- (x.out[x.out!=0])^exponent; x.out}


GetRconcavityValues <- function(pr.vec,n.s=100,r.val=-0.25,col=NULL) {
  pr.vec.copy <- pr.vec
  if (!any(! (unique(pr.vec) %in% round(seq(from=0,to=1,length=n.s+1),7)))) {
    pr.vec <- round(pr.vec * n.s) + 1
    min.prob <- min(pr.vec)
    pr.vec <- pr.vec - min(pr.vec)+1
    max.prob <- max(pr.vec)
    barplot(pos.exp(prop.table(tabulate(pr.vec,nbins=max.prob)),exponent=r.val),names.arg=as.character(seq(from=(min.prob-1)/n.s,to=(max.prob+min.prob-2)/n.s,length=max.prob)),col=col) # names(table(pr.vec.copy))
  } else {
    cat("ERROR","\n")
  }

}
# GetRconcavityValues(pr.vec=selection.probabilities[478,],r.val=-1)

# barplot(prop.table(table(x)))
# tabulate(pr.vec,nbins=n.s+1)


RunEmma <- function(gwas.obj,tr.n,cov.cols=integer(0),snp.name="",K.user=NULL) {
#gwas.obj=GWAS.obj;cov.cols=integer(0);snp.name=""

# OBJECTIVE :
# * REML estimates of genetic- and residual variance components, using the functions in emma.r
# *
# INPUT :
# * gwas.obj
# * tr.n : the trait number, i.e. column number in gwas.obj$pheno
# * cov.cols (optional) : column numbers (again in  gwas.obj$pheno) of covariates to be used
# * snp.name (optional) : name of a snp marker to be included as "covariate"
# * K.user (optional) : kinship matrix. Default is gwas.obj$kinship
# OUTPUT :
#  a vector with (in that order) the genetic- and residual variance
# TO DO : fixed effects p-values

  if (!is.null(K.user)) {gwas.obj$kinship <- K.user}

  if (dim(gwas.obj$kinship)[1] < gwas.obj$n) {

    subSet <- as.character(gwas.obj$pheno$genotype[!is.na(gwas.obj$pheno[,tr.n])])
  	K.emma <- gwas.obj$kinship[subSet,subSet]

 		rownames(K.emma) <- row.names(gwas.obj$pheno)[!is.na(gwas.obj$pheno[,tr.n])]
	  colnames(K.emma) <- row.names(gwas.obj$pheno)[!is.na(gwas.obj$pheno[,tr.n])]

	  temp.K <- K.emma

  } else {

    K.emma           <- gwas.obj$kinship[as.character(gwas.obj$pheno$genotype),as.character(gwas.obj$pheno$genotype)] #
	  rownames(K.emma) <- row.names(gwas.obj$pheno)
	  colnames(K.emma) <- row.names(gwas.obj$pheno)

  }

	n.var <- 1 + length(cov.cols)

  if (dim(gwas.obj$kinship)[1] < gwas.obj$n) {

  	if (n.var==1) {
  	  nonmissing <- which(!is.na(gwas.obj$pheno[,tr.n]))
  	} else {
  	  nonmissing <- which(apply(gwas.obj$pheno[!is.na(gwas.obj$pheno[,tr.n]),unique(c(tr.n,cov.cols))],
                                1,function(x){sum(!is.na(x))})==n.var)
  	}

	} else {

  	if (n.var==1) {
  	  nonmissing <- which(!is.na(gwas.obj$pheno[,tr.n]))
  	} else {
  	  nonmissing <- which(apply(gwas.obj$pheno[,unique(c(tr.n,cov.cols))],1,function(x){sum(!is.na(x))})==n.var)
  	}
  	temp.K <- K.emma[nonmissing,nonmissing]

	}

	temp.trait <- gwas.obj$pheno[nonmissing,tr.n]

	names(temp.trait) <- colnames(temp.K)

	n.ind  <- length(nonmissing)

	Xo   <- rep(1,n.ind)

	if (n.var>1) {
	  Xo <-cbind(Xo,gwas.obj$pheno[nonmissing,cov.cols])
	}
  if (snp.name!="") {
	  Xo <-cbind(Xo,as.numeric(gwas.obj$markers[snp.name,gwas.obj$pheno$genotype][nonmissing]))
	}

	Xo <- as.matrix(Xo)

	reml.obj <- emma.REMLE(temp.trait,Xo,temp.K)

return(list(varcomp=c(reml.obj$vg,reml.obj$ve),temp.K=temp.K))
}
# tests:
# RunEmma(gwas.obj=GWAS.obj,9)
# RunEmma(gwas.obj=GWAS.obj,9,10:11)
# RunEmma(gwas.obj=GWAS.obj,9,10:11,"m100")

RunEmmaWithWeights <- function(gwas.obj,tr.n,emma.weights=NULL,cov.cols=integer(0),snp.name="") {
#gwas.obj=GWAS.obj;cov.cols=integer(0);snp.name=""

# OBJECTIVE :
# * REML estimates of genetic- and residual variance components, using the functions in emma.r
# *
# INPUT :
# * gwas.obj
# * tr.n : the trait number, i.e. column number in gwas.obj$pheno
# * cov.cols (optional) : column numbers (again in  gwas.obj$pheno) of covariates to be used
# * snp.name (optional) : name of a snp marker to be included as "covariate"
# * emma.weights : if the phenotypic observations have different variances, and this is to be taken into account, emma.weights should be the vector
#                    of inverses of the standard deviations (not variances!)
# OUTPUT :
#  a vector with (in that order) the genetic- and residual variance
# TO DO : fixed effects p-values

  if (dim(gwas.obj$kinship)[1] < gwas.obj$n) {
  	K.emma      <- gwas.obj$kinship[gwas.obj$pheno$genotype[!is.na(gwas.obj$pheno[,tr.n])],gwas.obj$pheno$genotype[!is.na(gwas.obj$pheno[,tr.n])]]
 		rownames(K.emma) <- row.names(gwas.obj$pheno)[!is.na(gwas.obj$pheno[,tr.n])]
	  colnames(K.emma) <- row.names(gwas.obj$pheno)[!is.na(gwas.obj$pheno[,tr.n])]
	  temp.K <- K.emma  # gwas.obj$K[nonmissing,nonmissing]
  } else {
  	K.emma      <- gwas.obj$kinship[gwas.obj$pheno$genotype,gwas.obj$pheno$genotype]
	  rownames(K.emma) <- row.names(gwas.obj$pheno)
	  colnames(K.emma) <- row.names(gwas.obj$pheno)
  }
	n.var <- 1 + length(cov.cols)

  if (dim(gwas.obj$kinship)[1] < gwas.obj$n) {
  	if (n.var==1) {
  	  #nonmissing <- which(!is.na(gwas.obj$pheno[!is.na(gwas.obj$pheno[,tr.n]),tr.n]))
  	  nonmissing <- which(!is.na(gwas.obj$pheno[,tr.n]))
  	} else {
  	  nonmissing <- which(apply(gwas.obj$pheno[!is.na(gwas.obj$pheno[,tr.n]),unique(c(tr.n,cov.cols))],1,function(x){sum(!is.na(x))})==n.var)
  	}
	} else {
  	if (n.var==1) {
  	  nonmissing <- which(!is.na(gwas.obj$pheno[,tr.n]))
  	} else {
  	  nonmissing <- which(apply(gwas.obj$pheno[,unique(c(tr.n,cov.cols))],1,function(x){sum(!is.na(x))})==n.var)
  	}
  	temp.K <- K.emma[nonmissing,nonmissing]  # gwas.obj$K[nonmissing,nonmissing]
	}
	temp.trait <- gwas.obj$pheno[nonmissing,tr.n]
	#temp.K <- K.emma[nonmissing,nonmissing]  # gwas.obj$K[nonmissing,nonmissing]
	names(temp.trait) <- colnames(temp.K)
	n.ind  <- length(nonmissing)

  if (is.null(emma.weights)) {
	  Xo   <- rep(1,n.ind)
  } else {
	  Xo   <- emma.weights
  }
	if (n.var>1) { # also multiply the covariates with the weights
	  #Xo <-cbind(Xo,gwas.obj$pheno[nonmissing,cov.cols])
	  Xo <-cbind(Xo,diag(emma.weights) %*% as.matrix(gwas.obj$pheno[nonmissing,cov.cols]))
	}
  if (snp.name!="") { # also multiply the marker-scores with the weights
	  #Xo <-cbind(Xo,as.numeric(gwas.obj$markers[snp.name,gwas.obj$pheno$genotype][nonmissing]))
	  Xo <-cbind(Xo,emma.weights * as.numeric(gwas.obj$markers[snp.name,gwas.obj$pheno$genotype][nonmissing]))
	}

	Xo <- as.matrix(Xo)

	reml.obj <- emma.REMLE(temp.trait,Xo,temp.K)

return(c(reml.obj$vg,reml.obj$ve))
}


PlinkToGwasObj <- function(fam.file,map.file,plink.name) {
# OBJECTIVE :
# build a gwas.obj, from a plink .map and .fam file. The plink binaries must also exist
# (under the names plink.name.bed etc)
#
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *
  stopifnot (file.exists(paste(plink.name,".bim",sep="")))
  stopifnot (file.exists(paste(plink.name,".bed",sep="")))
  stopifnot (file.exists(paste(plink.name,".fam",sep="")))

  map <- read.table(map.file)
  fam <- read.table(fam.file)

  pheno <- data.frame(genotype=fam[,2],family=fam[,1],sex=fam[,5],disease=fam[,6])
  pheno$genotype <-as.character(pheno$genotype)
  pheno$family   <-as.character(pheno$family)
  pheno$disease  <- pheno$disease - 1
  mapframe <- data.frame(chromosome=map[,1],position=map[,4],snp.name=as.character(map[,2]))
  mapframe$snp.name <- as.character(mapframe$snp.name)
  row.names(mapframe) <- mapframe$snp.name

  nchr            <- length(unique(mapframe[,1]))
  chr.lengths     <- as.integer(table(mapframe$chromosome))
  #chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])


  if (nchr > 1) {
    chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])
    chr.lengths.bp  <- c(0,mapframe$position[cumsum(chr.lengths)[1:(nchr-1)]])
    cumposition     <- mapframe$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
  } else {
    cumposition     <- mapframe$position
    chr.pos         <- 0
    #chr.lengths.bp  <- nrow(marker.object)
  }
  mapframe <- cbind(mapframe,cum.position=cumposition)

  chromosomes <- sort(unique(mapframe$chromosome))

  gwas.obj <- list(map=mapframe,markers=data.frame(),pheno=pheno,external=list(plink.name=plink.name),plant.names=pheno$genotype,chromosomes=chromosomes,
                   nchr=length(chromosomes),n=nrow(fam),N=nrow(map),pca=data.frame(pca1=rep(0,nrow(fam)))) # chr.lengths

return(gwas.obj)
}

PlinkAssoc <- function(gwas.obj,pheno.file="",output.file="plinkoutput") {
# OBJECTIVE :
# association mapping using plink (armitage trend test or linear regression)
#
# INPUT :
# * gwas.obj : object created with the function PlinkToGwasObj
# * name of the output file (without extension)
# OUTPUT : none

  #stopifnot (exists(gwas.obj$external$plink.name))
  #stopifnot (file.exists(pheno.file))
  stopifnot (file.exists(paste(gwas.obj$external$plink.name,".bim",sep="")))
  stopifnot (file.exists(paste(gwas.obj$external$plink.name,".bed",sep="")))
  stopifnot (file.exists(paste(gwas.obj$external$plink.name,".fam",sep="")))
    # TO DO : WINDOWS/LINUX
  #stopifnot (file.exists("plink.exe"))
  if (pheno.file=="") {
    command <- paste("plink --bfile",gwas.obj$external$plink.name,"--assoc --out",output.file)
  } else {
    stopifnot (file.exists(pheno.file))
    command <- paste("plink --bfile",gwas.obj$external$plink.name,"--pheno",pheno.file,"--assoc --out",output.file)
  }
  cat(command,"\n")
  system(command)
}

ReadPlinkAssocFile <- function(assoc.file) {
  a <- read.table(file=assoc.file,header=T)
  a <- a[,c(2,7,4:6,8,9)]
  names(a) <- c("marker","major_allele","minor_allele","minorfreq.cases","minorfreq.controls","stat","pvalue")
  a$pvalue[is.na(a$pvalue)] <- 1
return(a)
}

AddBinaryPlinkFilesToGwasObj <- function(gwas.obj,plink.name) {
# to do : option make.files = F ; when true call MakePlinkFiles
#
# OBJECTIVE : 'associate' a set of 3 binary plink-files (bim/bed/fam) with the GWAS-object gwas.obj
#              these files should reside in the current working directory
#
# INPUT : gwas.obj and plink.name, the name of the 3 files (without extension)
# OUTPUT : the same object gwas.obj, but with an extra field "plink.name" in gwas.obj$external
#
  gwas.obj$external$plink.name <- plink.name
  stopifnot (file.exists(paste(plink.name,".bim",sep="")))
  stopifnot (file.exists(paste(plink.name,".bed",sep="")))
  stopifnot (file.exists(paste(plink.name,".fam",sep="")))
return(gwas.obj)
}





FaST_LMM <- function(gwas.obj,trait.nr,cov.cols=integer(0),full.lm=F,sim.file=character(0),
                     linreg=F,excludebyposition=0,ibs.kinship=F,
                     cov.type="RRM",kinship.snp.file=character(0),kinship.file=character(0),
                     simOutFile=character(0),output.file="output",show.call=F,
                     prefix="") {
#
# OBJECTIVE :
# *
# *
# INPUT :
# *
# * kinship.snp.file : file (including extension!) with the names of the snps  that are used for calculation of the kinship matrix
# * sim.file : file containing the kinship matrix to be used
#
# * cov.type :  "RRM" or "COVARIANCE" (ignored if sim.file is not character(0))
# *
# *
# OUTPUT : none; in contrast to the scan_GLS function, we only execute the program; reading the results (using the function ReadFastLmmResult) happens in the main script
# *
# *
# *
#

  # if (ibs.kinship) {}
  # only if .......
  #WriteFastLmmKinship(gwas.obj$kinship,tr.n=tr.n,file="temp_kinship.txt")
  # to do  : take into account covariates ....
  #command <- paste(command,"-sim temp_kinship.txt")
  if (ibs.kinship) {
    WriteFastLmmKinship(gwas.obj$kinship,tr.n=tr.n,file="temp_kinship.txt")
    command <- paste(command,"-sim temp_kinship.txt")
  }

  # to do  : input/output kinship matrix  (sim.file)
  #
  if (!(cov.type %in% c("RRM","COVARIANCE"))) {cov.type <- "RRM"}
  #stopifnot (exists(paste("plink.name",".bim",sep="")))
  #stopifnot(exists(gwas.obj$external$plink.name))
  stopifnot(length(gwas.obj$external$plink.name)==1)
  # to do : linux / windows specific:
  #stopifnot(file.exists("FastLmmC.exe"))
  #stopifnot(file.exists("libiomp5md.dll"))
  pheno.file <- names(gwas.obj$pheno)[trait.nr]
  MakeTfamFile(gwas.obj=gwas.obj,trait.nr  = trait.nr,file.name=pheno.file,file.type=1)

  # test if the plink file exists
  # test if the fast_lmm executable is in getwd()
  # To do: covariates ; loop over trait.nr 's

  #command.string      <- paste("fastlmmc",gwas.obj$external$bin.name,input.pheno,gwas.obj$external$kinship.name,varcomp.file,output.file,cov.string)

  os <- Sys.info()[['sysname']]
  if (!(os %in% c("Windows","Linux"))) {stop("The current version of the scripts require Windows or Linux.")}

  if (os=="Windows") {
    command.string      <- paste0(prefix,"FastLmmC -REML -ML -bfile ",gwas.obj$external$plink.name," -bfilesim ",
                                 gwas.obj$external$plink.name," -pheno ",paste0(pheno.file,".txt"),
                                 " -out ",output.file)
  }
  if (os=="Linux") {
    command.string      <- paste0(prefix,"fastlmmc -REML -ML -bfile ",gwas.obj$external$plink.name," -bfilesim ",
                                 gwas.obj$external$plink.name," -pheno ",paste0(pheno.file,".txt"),
                                 " -out ",output.file)
  }


  if (length(cov.cols)>0) {
    covar.file <- paste(names(gwas.obj$pheno)[trait.nr],"_covar",sep="")
    MakeTfamFile(gwas.obj,trait.nr  = cov.cols,file.name=covar.file,file.type=1,delimiter="\t")
    command.string      <- paste(command.string,"-covar",paste(covar.file,".txt",sep=""))
  }
  if (linreg) {
    sim.file=character(0);excludebyposition=0;extractSimTopK=character(0)
    command.string      <- paste(command.string,"-linreg")
  } else {
    if (cov.type == "COVARIANCE") {command.string      <- paste(command.string,"-simType COVARIANCE")}
    if (length(simOutFile)>0) {command.string      <- paste(command.string,"-simOut",simOutFile)}
    if (excludebyposition>0) {    command.string      <- paste(command.string,"-excludebyposition",as.character(excludebyposition))}
    if (length(kinship.file)>0) {command.string      <- paste(command.string,"-sim",kinship.file)}
    if (length(kinship.snp.file)>0) {
      m.snp <- nrow(read.table(kinship.snp.file))
      command.string      <- paste(command.string,"-extractSimTopK",kinship.snp.file,as.character(m.snp))
    }
    if (full.lm) {command.string   <- paste(command.string,"-simLearnType Full")}
  }
  if (show.call) {cat(command.string,"\n")}
cat(command.string,"\n")
system(command.string)
#cat(command.string,file="temp.bat")
#shell(cmd="temp.bat")
#shell(cmd="FastLMMc")
}

# FaST_LMM(gwas.obj= GWAS.obj,trait.nr=44,linreg=F,excludebyposition=0,kinship.snp.file=character(0),output.file="test_function")

ReplaceNaByOne <- function(x) {x[(is.na(x) | is.nan(x)) | x==-1]<-1; return(x)}
#testK <- FindOptimalM(gwas.obj=GWAS.obj,kinship.snp.file="best_snps_file.txt",plink.geno.name=GWAS.obj$external$plink.name,tr.n=9,m.start=100,m.step=100,inflation.threshold=1.2,ibs.kinship=T)
#gwas.obj=GWAS.obj;kinship.snp.file="best_snps_file.txt";plink.geno.name=GWAS.obj$external$plink.name;tr.n=9;m.start=100;m.step=100;inflation.threshold=1.2

ReadFastLmmKinship <- function(file.name="optimalkinship.txt",decimals=7) {
  M <- read.table(file.name,sep="\t",header=T) # row.names=1
return(round(M[,-1],decimals))
}

ReplaceStringInFile <- function(file.name,string.to.be.replaced,replacement="") {
# OBJECTIVE : replace a string in a text file by another string
# INPUT :
# OUTPUT : nothing returned; the text-file is rewritten
  con <- file(file.name)
  temp.string <- readLines(con = con, n=10^8)
  close(con)
  temp.string <- gsub(pattern=string.to.be.replaced, replacement=replacement, x=temp.string)
  cat(temp.string , file = file.name, sep = "\n")
}

ReadFastLmmResult <- function(file.name) {
  ReplaceStringInFile(file.name,"#IND")
  #M <- read.table(file.name,sep=",",header=T) # row.names=1
  Ncol <- count.fields(file=file.name,sep=',')[1]
  M <- read.table(file.name,sep=",",header=T,colClasses=c('character',rep('numeric',Ncol-1)))
  #M$WaldStat <- as.numeric(as.character(M$WaldStat))
  M$SNP <- as.character(M$SNP)
  M <- M[order(M$Chromosome,M$Position),]
  names(M)[1:5] <- c("snp","chromosome","GeneticDistance","position","pvalue")
  if ("WaldStat" %in% names(M)) {
    M$WaldStat <- (M$WaldStat)^2 / M$N
    #names(M)[which(names(M)=="GeneticDistance")] <- "stat"
  }
  # old code:
  #if ("WaldStat" %in% names(M)) {
  #  M$GeneticDistance <- (M$WaldStat)^2 / M$N
  #  names(M)[which(names(M)=="GeneticDistance")] <- "stat"
  #}
return(M)
}
#a <- ReadFastLmmResult("results/FT.output3_8.csv")


WriteFastLmmKinship <- function(gwas.obj,snp.subset=1:gwas.obj$N,file.name= "kinship.txt",tr.n=1) {
# N.B. assumes that the option accessions.as.families in the MakeTfamFile function was set to TRUE
# N.B. This matrix corresponds to only the nonmissing values; therefore, indicate the trait.number : tr.n
#       the complete matrix can always be obtained by choosing the genotype column (number 1, the default) as "trait"
  non.missing <- which(!is.na(gwas.obj$pheno[,tr.n]))
  n.ind <- length(non.missing) #nrow(gwas.obj$pheno)
  #ind.names <- as.character(non.missing) # as.character(1:n.ind)
  ind.names <- row.names(gwas.obj$pheno)[non.missing]
  ##cat(c("var",rep(ind.names,each=round(n.ind/gwas.obj$n))),file=file.name,sep="\t")
  mheader <- c("var",rep("",n.ind))
  mheader[-c(1)] <- paste(as.character(gwas.obj$pheno$genotype[non.missing]),ind.names)
  cat(mheader,file=file.name,sep="\t")
  cat("\n",file=file.name,append=T)
  if (length(snp.subset)!=gwas.obj$N) {
    temp.kinship <- IBS(gwas.obj$markers[snp.subset,])
    colnames(temp.kinship) <- gwas.obj$plant.names
    rownames(temp.kinship) <- gwas.obj$plant.names
    #a <- gwas.obj$kinship[gwas.obj$pheno$genotype,gwas.obj$pheno$genotype][non.missing,non.missing]
    temp.kinship <- temp.kinship[gwas.obj$pheno$genotype,gwas.obj$pheno$genotype][non.missing,non.missing]
  } else {
    temp.kinship <- gwas.obj$kinship[gwas.obj$pheno$genotype[non.missing],gwas.obj$pheno$genotype[non.missing]]
  }
  write.table(temp.kinship,sep="\t",row.names=mheader[-c(1)],col.names=F,quote=F,file=file.name,append=T)
}
#WriteFastLmmKinship(GWAS.obj,1:10,tr.n=tr.n)


MakeScanGlsKinship <- function(Sa2=1,Se2=1,K.matrix,plant.names,pheno.frame,tr.n=2) {
  # the following line is new : (22-2-2013)
  pheno.frame$genotype    <- as.character(pheno.frame$genotype)
	non.missing.pheno       <- !is.na(pheno.frame[,tr.n])
	n.ind                   <- sum(non.missing.pheno)
	geno.col.temp           <- pheno.frame$genotype[non.missing.pheno]
	TAU                     <- Sa2 * K.matrix[geno.col.temp,geno.col.temp] + Se2 * diag(n.ind)
	return(TAU)
}

# generalized version, with an error correlation matrix D.matrix
# D.matrix must be in the same format as K.matrix, i.e. of class matrix, and having row and column names (the genotype names)
MakeScanGlsKinship2 <- function(Sa2=1,Se2=1,K.matrix,D.matrix,plant.names,pheno.frame,tr.n=2) {
  pheno.frame$genotype    <- as.character(pheno.frame$genotype)
	non.missing.pheno       <- !is.na(pheno.frame[,tr.n])
	n.ind                   <- sum(non.missing.pheno)
	geno.col.temp           <- pheno.frame$genotype[non.missing.pheno]
	TAU                     <- Sa2 * K.matrix[geno.col.temp,geno.col.temp] + Se2 * D.matrix[geno.col.temp,geno.col.temp]
	return(TAU)
}


MakeGctaKinship <- function(Sa2.vector=1,Se2=1,gwas.obj,tr.n=2,asreml.components) {
	non.missing.pheno       <- !is.na(gwas.obj$pheno[,tr.n])
	n.ind                   <- sum(non.missing.pheno)
	geno.col.temp           <- gwas.obj$pheno$genotype[non.missing.pheno]

	TAU         <- diag(n.ind)
	for (chr in 1:length(Sa2.vector)) {
    TAU <- TAU + Sa2.vector[chr] * gwas.obj[[asreml.components[[chr]]]][geno.col.temp,geno.col.temp]#list.of.K.matrices[[chr]]
	}
	return(TAU)
}
#test.matrix2 <- MakeVcovMatrix2(varcomp.values[1,],varcomp.values[2,],GWAS.obj$kinship,plant.names=GWAS.obj$plant.names,GWAS.obj$pheno,9)
#write.table(test.matrix,file="vcov_test2.csv",sep=",",quote=F,col.names=F,row.names=F)

# an older, less efficient version :
#MakeVcovMatrix <- function(Sa2,Se2,K.matrix,plant.names,pheno.frame,tr.n=2) {
#		non.missing.pheno       <- !is.na(pheno.frame[,tr.n])
#		n.ind                   <- sum(non.missing.pheno)
#		geno.col.temp           <- pheno.frame$genotype[non.missing.pheno]
#
#        number.of.nonmissing <- aggregate(pheno.frame[,tr.n],by=list(ordered(pheno.frame$genotype)),FUN=function(x){sum(!is.na(x))})[match(plant.names,sort(plant.names)),2]
#		block.matrix<- matrix(0,n.ind,n.ind)
#		for (i in 1:n.ind) {
#		  block.matrix[i,] <- rep(as.numeric(GWAS.obj$kinship[which((pheno.frame$genotype[non.missing.pheno])[i]==plant.names)[1],]),times=number.of.nonmissing)
#		}
#		TAU         <- Sa2 * block.matrix + Se2 * diag(n.ind)
#		return(TAU)
#	}


FindOptimalM <- function(gwas.obj,kinship.snp.file,plink.geno.name,tr.n,m.start=100,m.step=100,sim.out=F,sim.out.name="optimalkinship.txt",
                         inflation.threshold=1.1,ibs.kinship=F,cov.cols=integer(0),
                         stop.at.local.minimum=F,only.stop.below.inflation.threshold=T,excludebyposition=0,group.snps=F) {   # ,excludebyposition=0
                        #,stop.below.GC.threshold=T

# OBJECTIVE :
# *
# *
# INPUT :
# * sim.out : write the new kinship matrix (GRM, based on the K selected snps) to a file ?  If TRUE, it is written to a file with name defined by sim.out.name
# * kinship.snp.file : name of the file to which the snps selected for the kinship matrix will be written to
# * stop.at.local.minimum : if TRUE, the search stops at the first step where the inflation-factor is higher than in the preceding step
# * inflation.threshold : if  stop.at.local.minimum  is FALSE (the default), the search stops once the inflation-factor is below inflation.threshold
# * only.stop.below.inflation.threshold  : when stop.at.local.minimum is TRUE, putting only.stop.below.inflation.threshold also to TRUE , will require that
#   the local minimum only 'counts' if at the same time the inflation-factor is below inflation.threshold

# OUTPUT :
# *
# *
# *


# TO DO : cov.cols unequal to integer(0)
#         make use of the FastLmm R-function

# input :

# GWAS.obj
#kinship.snp.file <- "best_snps_file.txt"
#tr.n <- 42
#plink.geno.name <- "LFN349acc_001"
#plink.geno.name <- "atwellData000"
#m.start <- 100
#m.step  <- 100
#inflation.threshold <- 1.15
# * sim.out : if TRUE, let fastlmm write the kinship matrix it has computed to a file
######

#if (sum(c(stop.at.local.minimum,stop.below.GC.threshold))==0) {stop.below.GC.threshold=T}

trait     <- names(gwas.obj$pheno)[tr.n]
m         <- m.start
m.vector  <- m
GC.vector <- numeric(0)
log.likelihood <- numeric(0)
stop.crit      <- F
old.inflation.factor <- 1000
cov.cols       <- integer(0)

# first, we do linear regression, without kinship
linreg.name <- paste("output/",trait,"_linreg.txt",sep="")
MakeTfamFile(gwas.obj,trait.nr  = tr.n,file.name=trait,file.type=1)
command <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",linreg.name,"-linreg")
cat(command,"\n")

# read the results of linear regression
system(command)
fastlmm.result <- read.table(file=linreg.name,header=T,sep="\t")
names(fastlmm.result)[1:5] <- c("snp","chromosome","GeneticDistance","position","pvalue")
fastlmm.result$snp <- as.character(fastlmm.result$snp)

if (group.snps) {
  fastlmm.result <- fastlmm.result[order(fastlmm.result$chromosome,fastlmm.result$position),]
}
row.names(fastlmm.result) <- fastlmm.result$snp

#ReadFastLmmResult(linreg.name)

########

if (group.snps) {
  # If group.snps=T, how should the grouping be done ?
  group.type  <- 1
  # 1- equally large groups of size group.size
  # 2- constructed with mapLD (to be done)
  # 3- based on genes (to be done)
  group.size  <- 30

  group.vector <- numeric(0)

  # INPUT : group.type, group.size
  # OUTPUT : K (number of groups)
  #          group.matrix ()
  #          group.list   ()
  #          group.vector ()
  #
  if (group.type==1) {
    group.list <- DefineGenomicRegions(gwas.obj$map$chromosome,block.size=group.size)
    group.matrix <- matrix(0,0,2)
    for (chr in 1:gwas.obj$nchr) {
      group.matrix <- rbind(group.matrix,group.list[[chr]]) #
    }
    K          <- nrow(group.matrix)
    if (length(unique(group.matrix[,2]-group.matrix[,1]))==1) {
      group.list <- data.frame(apply(group.matrix,1,function(x){(x[1]):(x[2])}))
    } else {
      group.list <- apply(group.matrix,1,function(x){(x[1]):(x[2])})
    }
    group.vector <- rep(0,gwas.obj$N)
    for (i in 1:K) {group.vector[group.matrix[i,1]:group.matrix[i,2]] <- i}
  }

  new.results <- matrix(-log10(fastlmm.result$pvalue))
  asd <- matrix(unlist(lapply(group.list,FUN=function(x){apply(matrix(new.results[x]),2,which.max)})),byrow=T,nrow=K)
  bestsnps <- sapply(1:length(asd),function(x){group.list[[x]][asd[x]]})
  bestsnps.frame <- data.frame(snp=as.character(fastlmm.result$snp[bestsnps]),pvalue=fastlmm.result$pvalue[bestsnps])
  bestsnps.frame <- bestsnps.frame[order(bestsnps.frame$pvalue),]
  bestsnps <-   as.character(bestsnps.frame$snp)
  #row.names(GWAS.obj$markers)
}
########

# define linreg.results as the list of best snps, of which we will take the K best
if (group.snps) {
  linreg.results <- bestsnps # fastlmm.result$snp[bestsnps]
} else {
  linreg.results <- fastlmm.result$snp[1:20000]
}

write.table(matrix(as.character(linreg.results[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)

while (!stop.crit) {

	final.results.name <- paste("output/",trait,"_final.txt",sep="")
  lik.results.name   <- paste("output/",trait,"_lik.txt",sep="")
	command            <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",final.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")
	command.lik        <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",lik.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")

  if (excludebyposition>0) {    command      <- paste(command,"-excludebyposition",as.character(excludebyposition))}

	if (ibs.kinship) {
	  WriteFastLmmKinship(gwas.obj,snp.subset=sort(which(row.names(gwas.obj$markers) %in% as.character(linreg.results[1:m]))),tr.n=tr.n,file="temp_kinship.txt")
	  command     <- paste(command,"-sim temp_kinship.txt")
	  command.lik <- paste(command.lik,"-sim temp_kinship.txt")
	}

	if (sim.out) {
	  command <- paste(command,"-simOut",sim.out.name)
	}

  cat(command,"\n")
  cat(command.lik,"\n")
	system(command)
	system(command.lik)

  ReplaceStringInFile(final.results.name,"#IND")
  ReplaceStringInFile(lik.results.name,"#IND")
	GWA.result <- read.table(file=final.results.name,header=T,sep="\t")
	GWA.result.lik <- read.table(file=lik.results.name,header=T,sep="\t")

	names(GWA.result)[1:5]  <- c("marker","chromosome","Position","position","pvalue")
	GWA.result              <- GWA.result[match(row.names(gwas.obj$markers),GWA.result$marker),]
	GC                      <- GenomicControlPvalues(pvals=GWA.result$pvalue,n.obs=sum(!is.na(gwas.obj$pheno[,tr.n])),n.cov=length(cov.cols))
	inflation.factor        <- GC[[2]]
  GC.vector               <- c(GC.vector,inflation.factor)
  log.likelihood          <- c(log.likelihood,GWA.result.lik$NullLogLike[1])

  # old version :
  #if (stop.at.local.minimum) {
	#  if (old.inflation.factor < inflation.factor | inflation.factor < inflation.threshold) {stop.crit <- T}
  #} else {
  #  if (inflation.factor < inflation.threshold) {stop.crit <- T}
  #}

  if (stop.at.local.minimum & only.stop.below.inflation.threshold) {
	  if (old.inflation.factor < inflation.factor & old.inflation.factor < inflation.threshold) {

      stop.crit <- T

      # extra: go one step back
      m <- m - m.step
      GC.vector <- GC.vector[1:(length(GC.vector)-1)]
      m.vector <- m.vector[1:(length(m.vector)-1)]
      log.likelihood <- log.likelihood[1:(length(log.likelihood)-1)]
      write.table(matrix(as.character(linreg.results[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)

	    command <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",final.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")
      if (excludebyposition>0) {    command      <- paste(command,"-excludebyposition",as.character(excludebyposition))}

    	if (ibs.kinship) {
    	  WriteFastLmmKinship(gwas.obj,snp.subset=sort(which(row.names(gwas.obj$markers) %in% as.character(linreg.results[1:m]))),tr.n=tr.n,file="temp_kinship.txt")
    	  command     <- paste(command,"-sim temp_kinship.txt")
    	}
    	if (sim.out) {
    	  command <- paste(command,"-simOut",sim.out.name)
    	}

      cat(command,"\n")
    	system(command)
      ReplaceStringInFile(final.results.name,"#IND")
    	GWA.result <- read.table(file=final.results.name,header=T,sep="\t")
    	names(GWA.result)[1:5]  <- c("marker","chromosome","Position","position","pvalue")
    	GWA.result              <- GWA.result[match(row.names(gwas.obj$markers),GWA.result$marker),]


    }
  }
  if (stop.at.local.minimum & !only.stop.below.inflation.threshold) {
	  if (old.inflation.factor < inflation.factor) {
      stop.crit <- T

      # extra: go one step back
      m <- m - m.step
      GC.vector <- GC.vector[1:(length(GC.vector)-1)]
      m.vector <- m.vector[1:(length(m.vector)-1)]
      log.likelihood <- log.likelihood[1:(length(log.likelihood)-1)]
      write.table(matrix(as.character(linreg.results[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)

	    command <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",final.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")
      if (excludebyposition>0) {    command      <- paste(command,"-excludebyposition",as.character(excludebyposition))}

    	if (ibs.kinship) {
    	  WriteFastLmmKinship(gwas.obj,snp.subset=sort(which(row.names(gwas.obj$markers) %in% as.character(linreg.results[1:m]))),tr.n=tr.n,file="temp_kinship.txt")
    	  command     <- paste(command,"-sim temp_kinship.txt")
    	}
    	if (sim.out) {
    	  command <- paste(command,"-simOut",sim.out.name)
    	}

      cat(command,"\n")
    	system(command)
      ReplaceStringInFile(final.results.name,"#IND")
    	GWA.result <- read.table(file=final.results.name,header=T,sep="\t")
    	names(GWA.result)[1:5]  <- c("marker","chromosome","Position","position","pvalue")
    	GWA.result              <- GWA.result[match(row.names(gwas.obj$markers),GWA.result$marker),]


    }
  }
  if (!stop.at.local.minimum) {
    if (inflation.factor < inflation.threshold) {stop.crit <- T}
  }

	if (!stop.crit) {
    m <- m + m.step
	  m.vector <- c(m.vector,m)
	  write.table(matrix(as.character(linreg.results[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)
	  old.inflation.factor <- inflation.factor
    cat("\n","\n","CURRENT INFLATION FACTOR =",inflation.factor,"\n","\n")
  }
} # end while


kinship.snps <-   read.table(file=kinship.snp.file,colClasses="character")

return(list(kinship.snps=kinship.snps,kinship.snp.file=kinship.snp.file,m.optimal=m,GC.vector=GC.vector,m.vector=m.vector,log.lik=log.likelihood))
}


FindOptimalM.orthogonal <- function(gwas.obj,kinship.snp.file,plink.geno.name,tr.n,m.step=5,sim.out=F,sim.out.name="optimalkinship.txt",
                                    inflation.threshold=1.1,cov.cols=integer(0),stop.at.local.minimum=F,excludebyposition=25000) {
#browser()

# goal
# (optional) :  start with a base set, e.g. of 139 snps
# find snps one by one, using linear regression. At Each step we do linear regression conditional on
#                            the snps found during earlier steps, + possibly the base set which is included all the time


#gwas.obj=GWAS.obj;plink.geno.name="rik"
#tr.n=8;m.step=5
#sim.out=F;sim.out.name="optimalkinship.txt"
#inflation.threshold=1.1
#
#stop.at.local.minimum=F;excludebyposition=25000

####################################

trait<- names(gwas.obj$pheno)[tr.n]

#linreg.name <- paste("output/",trait,"_linreg.txt",sep="")

covar.vector <- character(0)
covar.file   <- "output/covar_file"

stop.crit <- F
#old.inflation.factor <- 1000
cov.cols <- integer(0)

output.file         <- paste("output/",trait,".","output_selectK",".csv",sep="")
kinship.snp.file    <- paste("output/",trait,".","snp_selection",".txt",sep="")

n.var          <- ncol(gwas.obj$pheno)

while (!stop.crit) {
  if (length(covar.vector) > 0) {
    gwas.obj$pheno <- cbind(gwas.obj$pheno,t(as.matrix(gwas.obj$markers[covar.vector,gwas.obj$pheno$genotype])))
    cov.cols <- (n.var+1):(ncol(gwas.obj$pheno))
  }
  FaST_LMM(gwas.obj=gwas.obj,trait.nr=tr.n,cov.cols=cov.cols,linreg=T,output.file=output.file,show.call=T)
  gwas.obj$pheno <- gwas.obj$pheno[,1:n.var]
  smallest.pvalue <- as.numeric(strsplit(readLines(con = output.file,n=2)[[2]],split=",")[[1]][5])
  if (smallest.pvalue  < 0.05 / gwas.obj$N) {
    covar.vector   <- c(covar.vector,strsplit(readLines(con = output.file,n=2)[[2]],split=",")[[1]][1])
  } else {
    stop.crit <- T
  }
}

marker.selection <- integer(0)
m <- length(covar.vector)
m.vector <- m
if (m>0) {
  stop.crit <- F
  GWA.result <- ReadFastLmmResult(output.file)
  marker.selection <- c(covar.vector,paste("m",SelectMarkers(GWA.result$pvalue,gwas.obj,25000,max.number=2000),sep=""))
  #write.table(matrix(as.character(covar.vector)),quote=F,file=kinship.snp.file,row.names=F,col.names=F)
  write.table(matrix(as.character(marker.selection[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)

  FaST_LMM(gwas.obj=gwas.obj,trait.nr=tr.n,excludebyposition=excludebyposition,output.file=output.file,kinship.snp.file=kinship.snp.file,show.call=T)
  GWA.result <- ReadFastLmmResult(output.file)
  #GWA.result$pvalue[minor.allele.frequencies < maf] <- NA
  smallest.pvalue <- min(GWA.result$pvalue,na.rm=T)
  GC <- GenomicControlPvalues(pvals=GWA.result$pvalue,n.obs=sum(!is.na(gwas.obj$pheno[,tr.n])),n.cov=length(cov.cols))
  inflation.factor  <- GC[[2]]
  GC.vector  <- inflation.factor
  log.likelihood <- mean(GWA.result$NullLogLike)
  if (smallest.pvalue  > 0.05 / gwas.obj$N) {stop.crit <- T}
  if (smallest.pvalue  <= 0.05 / gwas.obj$N & inflation.factor <= inflation.threshold) {
    stop.crit <- T
  }
} else {
  kinship.snps=integer(0)
  kinship.snp.file=character(0)
  m.optimal=m
  GC.vector=NULL
  m.vector=m
  log.likelihood <- NULL
}


while (!stop.crit) {
   m <- m + m.step
   m.vector <- c(m.vector,m)
  write.table(matrix(as.character(marker.selection[1:m])),quote=F,file=kinship.snp.file,row.names=F,col.names=F)
  old.inflation.factor <- inflation.factor
  cat("\n","\n","CURRENT INFLATION FACTOR =",inflation.factor,"\n","\n")


  FaST_LMM(gwas.obj=gwas.obj,trait.nr=tr.n,excludebyposition=excludebyposition,output.file=output.file,kinship.snp.file=kinship.snp.file,show.call=T,simOutFile=sim.out.name)
	#final.results.name <- paste("output/",trait,"_final.txt",sep="")
  #lik.results.name <- paste("output/",trait,"_lik.txt",sep="")
	#command <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",final.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")
	#command.lik <- paste("FastLmmC -bfile",plink.geno.name,"-bfilesim",plink.geno.name,"-pheno",paste(trait,".txt",sep=""),"-out",lik.results.name,"-extractSimTopK",kinship.snp.file,m," -REML")
  GWA.result <- ReadFastLmmResult(output.file)
	GC         <- GenomicControlPvalues(pvals=GWA.result$pvalue,n.obs=sum(!is.na(gwas.obj$pheno[,tr.n])),n.cov=length(cov.cols))
	inflation.factor  <- GC[[2]]
  GC.vector  <- c(GC.vector,inflation.factor)
  log.likelihood <- c(log.likelihood,mean(GWA.result$NullLogLike))

  if (stop.at.local.minimum) {
	  if (old.inflation.factor < inflation.factor | inflation.factor < inflation.threshold) {stop.crit <- T}
  } else {
    if (inflation.factor < inflation.threshold) {stop.crit <- T}
  }

}

#kinship.snps <-   read.table(file=kinship.snp.file,colClasses="character")
if (m>0) {
  kinship.snps <- as.character(marker.selection[1:m])
}

return(list(kinship.snps=kinship.snps,kinship.snp.file=kinship.snp.file,m.optimal=m,GC.vector=GC.vector,m.vector=m.vector,log.lik=log.likelihood))

}

JoinSubsamples <- function(s1,s2) {

  s1$sub.samples <- cbind(s1$sub.samples,s2$sub.samples)
  for (i in 1:nrow(s1$sub.samples)) {
    nonzero <- sort(s1$sub.samples[i,][s1$sub.samples[i,]!=0])
    s1$sub.samples[i,] <- 0
    s1$sub.samples[i,1:length(nonzero)] <- nonzero
  }
  s1$sub.samples <- s1$sub.samples[,which(apply(s1$sub.samples,2,sum)>0)]

  s1$sub.sample.size   <- s1$sub.sample.size + s2$sub.sample.size
  s1$ind.population    <- sort(union(s1$ind.population,s2$ind.population))
  s1$ind.sample.size   <- s1$ind.sample.size + s2$ind.sample.size
  s1$number.of.nonmissing  <- s1$number.of.nonmissing + s2$number.of.nonmissing
  s1$acc.population        <- sort(union(s1$acc.population ,s2$acc.population ))
  s1$acc.sample.size       <- s1$acc.sample.size + s2$acc.sample.size

return(s1)
}
# typical use :
#temp.pheno <- GWAS.obj$pheno
#temp.pheno$disease[temp.pheno$disease==0] <- NA
#s1   <- DrawSubsamples(pheno.dataframe=temp.pheno,trn=tr.n,plant.names=GWAS.obj$plant.names,ns=NS,subsample.acc=TRUE)
#temp.pheno <- GWAS.obj$pheno
#temp.pheno$disease[temp.pheno$disease==1] <- NA
#s2   <- DrawSubsamples(pheno.dataframe=temp.pheno,trn=tr.n,plant.names=GWAS.obj$plant.names,ns=NS,subsample.acc=TRUE)
#asd <- JoinSubsamples(s1,s2)

GetSNPsInRegion <- function(gwas.obj,snp.number,region.size=5000){
  return(setdiff(which((abs(gwas.obj$map$position[snp.number] - gwas.obj$map$position) <= region.size) & (gwas.obj$map$chromosome==gwas.obj$map$chromosome[snp.number])),snp.number))
}

SelectMarkers <- function(pvalues,gwas.obj,window.size,min.p=0.1,max.number=length(pvalues)) {
# OBJECTIVE : given a (physical) map and a vector of p-values, construct a vector of markers (more precisely, marker numbers) that are most significant, and that are not in LD.
#   First, the most significant marker is selected, and a region around this marker is exlcuded.
#   Subsequent markers are found similarly, by selecting the most significant markers outside the 'LD-regions' of the markers that have already been selected.
#   The process stops if either these LD-regions completely cover the complete genome, or if there no markers left with p-value smaller than a threshold (default min.p=0.1)
# INPUT :
#   gwas.obj : a GWAS-object, with gwas.obj$N markers
#   pvalues : a vector of pvalues fo length gwas.obj$N, obtained for gwas.obj
#   window.size : this number indicates how many base -pairs on both sides of a selected marker are excluded
# OUTPUT :
#   the marker numbers, i.e. element numbers in the vector pvalues
#
  stop.search <- F
  available.markers <- rep(T,gwas.obj$N)
  available.markers[is.na(pvalues)] <- F   # missing pvalues are not taken into account
  selected.markers <- integer(0)
  m <- 0
  while (!stop.search) {
    smallest.p <- min(pvalues[available.markers])
    best.marker <- which.min(pvalues[available.markers])[1]  # N.B if several markers have the same (minimal) pvalue, we choose the first one
    best.marker <- ((1:gwas.obj$N)[available.markers])[best.marker]
    selected.markers <- c(selected.markers,best.marker)
    available.markers[c(GetSNPsInRegion(gwas.obj,best.marker,region.size=window.size),best.marker)] <- F
    if (sum(available.markers)==0 | smallest.p > min.p) {stop.search <- T}
    m <- m + 1
    if (m >= max.number) {stop.search <- T}
  }
return(sort(selected.markers)) # output : the marker numbers, i.e. element numbers in the vector pvalues
}
#qw <- SelectMarkers(GWA.result$pvalue,GWAS.obj,400000)


SelectMarkersMin <- function(pvalues,gwas.obj,window.size,min.p=0.1) {
# OBJECTIVE : given a (physical) map and a vector of p-values, construct a vector of markers (more precisely, marker numbers) that are LEAST significant, and that are not in LD.
#   First, the least significant marker is selected, and a region around this marker is exlcuded.
#   Subsequent markers are found similarly, by selecting the least significant markers outside the 'LD-regions' of the markers that have already been selected.
#   The process stops if either these LD-regions completely cover the complete genome, or if there no markers left with p-value larger than a threshold (default min.p=0.1)
# INPUT :
#   gwas.obj : a GWAS-object, with gwas.obj$N markers
#   pvalues : a vector of pvalues fo length gwas.obj$N, obtained for gwas.obj
#   window.size : this number indicates how many base -pairs on both sides of a selected marker are excluded
# OUTPUT :
#   the marker numbers, i.e. element numbers in the vector pvalues
#
  stop.search <- F
  available.markers <- rep(T,gwas.obj$N)
  available.markers[is.na(pvalues)] <- F   # missing pvalues are not taken into account
  selected.markers <- integer(0)
  while (!stop.search) {
    largest.p <- max(pvalues[available.markers])
    best.marker <- which.max(pvalues[available.markers])[1]  # N.B if several markers have the same (minimal) pvalue, we choose the first one
    best.marker <- ((1:gwas.obj$N)[available.markers])[best.marker]
    selected.markers <- c(selected.markers,best.marker)
    available.markers[c(GetSNPsInRegion(gwas.obj,best.marker,region.size=window.size),best.marker)] <- F
    if (sum(available.markers)==0 | largest.p < min.p) {stop.search <- T}
  }
return(sort(selected.markers)) # output : the marker numbers, i.e. element numbers in the vector pvalues
}


ComputeEigenstratMatrix <- function(gwas.obj,maf=0,chr=0) {
  X <- as.matrix(gwas.obj$markers)
  if (chr!=0){
    X <- X[gwas.obj$map$chromosome==chr,]
  }
  minor.allele.frequencies <- apply(X,1,function(x){mean(as.numeric(x))})
  X <- X[minor.allele.frequencies >= maf & minor.allele.frequencies <= (1-maf),]
  minor.allele.frequencies <- minor.allele.frequencies[minor.allele.frequencies >= maf & minor.allele.frequencies <= (1-maf)]
  stdev <- sqrt(minor.allele.frequencies*(1-minor.allele.frequencies))
  X <- X - matrix(rep(minor.allele.frequencies,gwas.obj$n),ncol=gwas.obj$n)
  for (j in 1:gwas.obj$n) {X[,j] <- X[,j]/stdev}
  M <- t(X) %*% X / nrow(X)
return(M)
}
#cov.matrix <- ComputeEigenstratMatrix(GWAS.obj)
#write.table(cov.matrix,file="D:/willem/statistical_genetics_large_files/arabidopsis_data/LFNdata/eigenstrat_matrix_LFN.csv",quote=F,sep=",",row.names=FALSE,col.names=GWAS.obj$plant.names)

QQplotPvalues <- function(pvalues) {
# makes a qq-plot of observed versus expected LOD-scores
# The code was taken from Segura et al (2012)
    pvalues <- na.omit(pvalues)
		e<--log10(ppoints(length(pvalues)))
		o<--log10(sort(pvalues))
		maxp<-ceiling(max(o))
		plot(e,o,type='b',pch=20,cex=0.8,col=1,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),
         xlim=c(0,max(e)+1),ylim=c(0,maxp),main=paste('QQ-plot'))
		abline(0,1,col="dark grey")
}
#GWA.result$pvalue

QQplotPvalues2 <- function(pvalues,main.title='QQ-plot',file.name="") {
# makes a qq-plot of observed versus expected LOD-scores
# The code was taken from Segura et al (2012) and adapted
    pvalues <- na.omit(pvalues)
		e       <- -log10(ppoints(length(pvalues)))
		o       <- -log10(sort(pvalues))
		maxp    <- ceiling(max(o))

    if (file.name=="") {
		plot(e,o,type='b',pch=20,cex=0.9,col=1,xlab=expression(Expected~~-log[10](p)),
                                           ylab=expression(Observed~~-log[10](p)),
                                           xlim=c(0,max(e)+1),ylim=c(0,maxp),
                                           main=main.title)
		                                       abline(0,1,col="blue")
    } else {
      #jpeg(file=file.name,quality=100)
      png(file=file.name)
  		plot(e,o,type='b',pch=20,cex=0.9,col=1,xlab=expression(Expected~~-log[10](p)),
                                             ylab=expression(Observed~~-log[10](p)),
                                             xlim=c(0,max(e)+1),ylim=c(0,maxp),
                                             main=main.title)
  		                                       abline(0,1,col="blue")
      dev.off()
    }
}
#GWA.result$pvalue


# create_gwasObj_from_files
genstatToGwasObj <- function(geno.file,map.file,pheno.file,kin.file,description="mydata") {
# OBJECTIVE :
# build a gwas.obj, from text-files containing the marker-scores, map (chr + base pairs or centi-Morgan), phenotypic data and the kinship matrix.
# Check the example files (geno.file="sim_locfile.csv",map.file="mapfile.txt",pheno.file="pheno.csv",kin.file="kinship.csv") for the required format.
# The files can, but not have to be created using genstat

  stopifnot (file.exists(geno.file))
  stopifnot (file.exists(pheno.file))
  stopifnot (file.exists(map.file))
  stopifnot (file.exists(kin.file))

  map <- read.table(map.file,header=T)
  markers <- read.table(geno.file,sep=",",header=T,row.names=1)
  pheno <-  read.table(pheno.file,sep=",",header=T)
  kin    <- as.matrix(read.table(kin.file,sep=",",header=T))
  colnames(kin) <- names(markers)
  rownames(kin) <- names(markers)

  names(pheno)[1] <- "genotype"
  pheno$genotype <-as.character(pheno$genotype)

  mapframe <- data.frame(chromosome=map[,1],position=map[,2],snp.name=as.character(map[,3]))
  mapframe$snp.name <- as.character(mapframe$snp.name)
  row.names(mapframe) <- mapframe$snp.name

  nchr            <- length(unique(mapframe[,1]))
  chr.lengths     <- as.integer(table(mapframe$chromosome))

  if (nchr > 1) {
    chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])
    chr.lengths.bp  <- c(0,mapframe$position[cumsum(chr.lengths)[1:(nchr-1)]])
    cumposition     <- mapframe$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
  } else {
    cumposition     <- mapframe$position
    chr.pos         <- 0
    chr.lengths.bp  <- nrow(markers)
  }
  mapframe <- cbind(mapframe,cum.position=cumposition)

  chromosomes <- sort(unique(mapframe$chromosome))
  N <- nrow(markers)
  n <- ncol(markers)

  bin.file.name <- paste(description,".bin",sep="")
  kinship.name <- paste(description,".csv",sep="")
  MakeBin(kinship=kin,plant.names=names(markers),pheno=pheno,csv.file.name=geno.file,bin.file.name=bin.file.name)

  gwas.obj <- list(map=mapframe,markers=markers,pheno=pheno,kinship=kin,external=list(bin.name=bin.file.name,kinship.name=kinship.name),plant.names=names(markers),chromosomes=chromosomes,
                   nchr=length(chromosomes),n=n,N=N,kinship.asreml=MakeKinshipAsreml(kin),chr.lengths.bp=chr.lengths.bp,genes=data.frame())

return(gwas.obj)
}

createZmatrixForIncompleteBlockDesign <- function(rep.vector,block.vector) {
# no missing values
# rep.vector and block.vector must be of equal length
  if (!is.factor(rep.vector)) {rep.vector <- factor(rep.vector)}
  if (!is.factor(block.vector)) {block.vector <- factor(block.vector)}
  stopifnot(length(rep.vector)==length(block.vector) )
  stopifnot(sum(is.na(rep.vector))==0 & sum(is.na(block.vector))==0)
  stopifnot(1 %in% levels(rep.vector) & 1 %in% levels(block.vector))
  new.factor <- rep.vector:block.vector
  temp.data <- data.frame(y=1:length(rep.vector),new.factor=new.factor) # rep.vector=rep.vector,block.vector=block.vector
  Zmatrix <- model.matrix(as.formula(y ~ new.factor),data=temp.data)#(GWAS.obj$pheno[,2] : GWAS.obj$pheno[,3])
  Zmatrix[,1] <- 0
  Zmatrix[rep.vector==1 & block.vector==1,1] <- 1
return(list(Zmatrix=Zmatrix,nested.block.factor=new.factor))
}

createZmatrixForFactor <- function(f.vector,factor.name="",ordered.factor=F) {
#browser()
# no missing values!!
#f.vector=gwas.obj$pheno$genotype;f.levels=gwas.obj$plant.names;ordered.factor=T;factor.name=""
  if (sum(is.na(f.vector))>0) {
    stop("Missing values in factor not allowed.")
  }
  if (class(f.vector)=="data.frame") {f.vector <- f.vector[,ncol(f.vector)]}
  if (!is.factor(f.vector)) {
    f.levels <- unique(as.character(sort(f.vector)))
    f.vector <- factor(f.vector,levels=f.levels,ordered=ordered.factor)
  }
  temp.data <- data.frame(y=1:length(f.vector),new.factor=f.vector)
  Zmatrix <- model.matrix(as.formula(y ~ new.factor),data=temp.data)
  Zmatrix[,1] <- 0
  Zmatrix[f.vector==levels(f.vector)[1],1] <- 1
  colnames(Zmatrix) <- paste(factor.name,levels(f.vector),sep="")
return(Zmatrix)
}
#createZmatrixForFactor(GWAS.obj$pheno$ronde,factor.name="ronde")
# Z.geno2 <- createZmatrixForFactor(gwas.obj$pheno$genotype,f.levels=gwas.obj$plant.names,T)


WriteGenstatCsv <- function(X,file.name="genstatData.csv") {
  factor.cols <- which(lapply(X,FUN=class)=="factor")
  names(X)[factor.cols] <- paste(names(X)[factor.cols],"!",sep="")
  write.csv(X,file=file.name,na="*",row.names=F,quote=F)
}

GeneralizedHeritability <- function(reml.obj,K.matrix) {
# OBJECTIVE : calculate the generalized heritability as proposed by Oakey et al (2006, TAG)
# INPUT :
# * reml.obj : asreml output obtained after fitting the variance components with asreml
#   The first element in summary(reml.obj)$varcomp$component  must come from 'genotype'; the last one is the residual error
# * K.matrix : the (unscaled) genetic variance times the kinship matrix of the genotypes. So if K is the kinship matrix and g ~ N(0,\sigma_g^2 K)
#              then K.matrix should be \hat \sigma_g^2 K, where \hat \sigma_g^2 is an estimate of \sigma_g^2
# OUTPUT (H2) :  the generalized heritability as proposed by Oakey et al
#
# EXAMPLE :
#  K.matrix <- summary(reml.obj)$varcomp$component[1] * GWAS.obj$kinship
#  GeneralizedHeritability(reml.obj,K.matrix)
#
  reml.pred <- predict(reml.obj, classify="genotype",vcov=T)

  A <- reml.pred$predictions$vcov
  #evs <- eigen(diag(ncol(A)) - solve(K.matrix) %*% A )$values
  # a small imaginary part may appear, due to numerical inaccuracies ?
  evs <- Re(eigen(diag(ncol(A)) - solve(K.matrix) %*% A )$values)
  H2 <- mean(evs[evs>0])
  return(H2)
}

SimulatePhenoGivenGenoAlhpaDesignWithChangingPanel <- function(gwas.obj,pheno.list,n.sim=10,rel.effects=c(0.05,0.1),
                                                               sigmaG2,sigmaE2,sigmaB2=0,add.pca=F,n.pred=0,
                                                               n.corrupt.geno=0,n.left=1,use.asreml=T,no.BLUPs=F,line.heritability=F,
                                                               MSE.method=F) {  # ,sub.set=GWAS.obj$plant.names
# new version: qtl-locations are drawn again for each simulation
# simulation of a phenotype given real genotypic GWAS data.
# An (already generated) alpha design is assumed
# INPUT :
# * sigmaG2 = the (scaled) polygenic variance
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
#   columns 2 and 3 of gwas.obj$pheno must be replicate and block
# * n.sim : number of simulated data sets
# * different.names : when TRUE, the simulated traits will have col-names sim1, sim2,...
#                    otherwise they will all have col-name pheno
# * n.pred : number of unobserved genotypes that are simulated and predicted
# * n.corrupt.geno : for this number of genotypes, leave n.left replicates nonmissing
# * n.left :
# * use.asreml :
# * no.BLUPs : calculate BLUPs for the observed individuals ?
# * MSE.method : perform also mixed model based estimation of h2 based on the means, but fixing sigma_e^2 to MS(Res)
#
# OUTPUT :
# *
# *

# gwas.obj=GWAS.obj;n.sim=length(pheno.list);rel.effects=c(0.05,0.1); n.pred=n.predict; use.asreml=F

# gwas.obj=GWAS.obj; n.sim=length(pheno.list) ; n.pred=n.predict
  if (is.null(gwas.obj$kinship2)) {
    gwas.obj$kinship2 <- gwas.obj$kinship
    gwas.obj$kinship.asreml2 <- MakeKinshipAsreml(gwas.obj$kinship2)
  }

  pheno.list.out <- pheno.list

  n.qtl <- length(rel.effects)
  n.rep <- nlevels((pheno.list[[1]])$replicate) # N.B. assumes equal numbers of replicates in the different simulations!!

  mafs <- list()
  candidates <- list()
  candidates.accessions <- unique(unlist(lapply(pheno.list,function(x){x$genotype})))
  for (j in gwas.obj$chromosomes) {  # N.B. GWAS.obj$chromosomes should be labeled 1,2,3...
    mafs[[j]] <- apply(gwas.obj$markers[gwas.obj$map$chromosome==j,candidates.accessions],1,mean)
    candidates[[j]] <- ((1:gwas.obj$N)[gwas.obj$map$chromosome==j])[which(mafs[[j]] > .1 & mafs[[j]] < .9)]
  }

  pos.frame    <- as.data.frame(matrix(0,n.qtl,n.sim))
  names(pos.frame) <- paste("sim",as.character(1:n.sim),sep="")
  chr.frame    <- pos.frame
  effect.frame <- pos.frame

  Z.block      <- createZmatrixForIncompleteBlockDesign(rep.vector=(pheno.list[[1]])$replicate,block.vector=(pheno.list[[1]])$block)
  block.factor <- Z.block$nested.block.factor
  Z.block      <- Z.block$Zmatrix

  explained.phenotypic.variance <- rep(NA,n.sim)

  qtl.variance.exact  <- rep(NA,n.sim)
  qtl.variance.independent <- rep(NA,n.sim)

  mixed.model.h2.estimate.1stage <- rep(NA,n.sim)  # mixed model estimate of h2, based on the replicates; calculated with emma
  mixed.model.h2.estimate.1stageB<- rep(NA,n.sim)  # mixed model estimate of h2, based on the replicates; calculated with asreml
  mixed.model.h2.estimate.2stage <- rep(NA,n.sim)  # mixed model estimate of h2, based on the genotypic means; calculated with emma
  mixed.model.h2.estimate.2stageB<- rep(NA,n.sim)  # mixed model estimate of h2, based on the genotypic means; calculated with asreml
  #
  mixed.model.h2.estimate.2stageBextra<- rep(NA,n.sim)  # mixed model estimate of h2, based on the genotypic means + MS(Res); calculated with asreml
  mixed.model.sg.estimate.2stageBextra<- rep(NA,n.sim)  # mixed model estimate of sg, based on the genotypic means + MS(Res); calculated with asreml
  mixed.model.se.estimate.2stageBextra<- rep(NA,n.sim)  # mixed model estimate of se, based on the genotypic means + MS(Res); calculated with asreml
  if (MSE.method) {stopifnot(use.asreml)}

  anova.h2.estimate <- rep(NA,n.sim)  # the plot heritability , e.g. MS(G) / (MS(G) + MS(E))
  anova.h2.estimate2 <- rep(NA,n.sim) # the ANOVA estimate , calculated with the Heritability3 function

  sigmaG2.estimate.1 <- rep(NA,n.sim) # mixed model estimate of genetic variance, based on the replicates;
  sigmaG2.estimate.2 <- rep(NA,n.sim) # mixed model estimate of genetic variance, based on the genotypic means
  sigmaE2.estimate.1 <- rep(NA,n.sim)
  sigmaE2.estimate.2 <- rep(NA,n.sim)

  conf.int1.left.1stage <- rep(NA,n.sim)   # left bound of the conf. interval for the mixed model based estimate heritability; based on the replicates. Only if asreml is used
  conf.int1.right.1stage <- rep(NA,n.sim)  # right bound of the conf. interval for the mixed model based estimate heritability; based on the replicates. Only if asreml is used
  conf.int1.left.2stage <- rep(NA,n.sim)   # left bound of the conf. interval for the mixed model based estimate heritability; based on the genotypic means. Only if asreml is used
  conf.int1.right.2stage <- rep(NA,n.sim)  # right bound of the conf. interval for the mixed model based estimate heritability; based on the genotypic means. Only if asreml is used

  conf.int2.left.1stage <- rep(NA,n.sim)   # the same as above, but intervals back transformed after using the delta method and intervals  for log(\sigma_g^2 / \sigma_e^2)
  conf.int2.right.1stage <- rep(NA,n.sim)
  conf.int2.left.2stage <- rep(NA,n.sim)
  conf.int2.right.2stage <- rep(NA,n.sim)

  conf.int3.left.1stage <- rep(NA,n.sim)   # unused ?
  conf.int3.right.1stage <- rep(NA,n.sim)
  conf.int3.left.2stage <- rep(NA,n.sim)
  conf.int3.right.2stage <- rep(NA,n.sim)

  conf.int.left.broad.sense  <- rep(NA,n.sim) # confidence interval for the broad-sense ANOVA heritability; only if line.heritability = F
  conf.int.right.broad.sense <- rep(NA,n.sim)

  K.sd.off.diagonal <- rep(NA,n.sim)

  breeding.values     <- list()              # list of data-frames with the BLUPs for the genetic effects of the observed individuals (+ a column for the real effects ?)
  breeding.values.MSE1 <- rep(NA,n.sim)      # the corresponding MSE's for 1= 1-stage = for BLUPs based on the means; 2 = 2-stage = for BLUPs based on the replicates
  breeding.values.MSE2 <- rep(NA,n.sim)

  # genomic selection : prediction of genotypic effect for N.pred
  GS                   <- list()           # list of data-frames with the predicted genetic effects for the unobserved individuals
  GS.MSE1 <- rep(NA,n.sim)
  GS.MSE2 <- rep(NA,n.sim)
  average.pred.error1 <- rep(NA,n.sim)
  average.pred.error2 <- rep(NA,n.sim)

  for (trait in 1:n.sim) {

    cat('\n','simulation ',trait,'\n\n')

    available.genotypes <- unique((pheno.list[[trait]])$genotype)
    n.available         <- length(available.genotypes)

    K <- gwas.obj$kinship[available.genotypes,available.genotypes]
    K.transform <- KinshipTransform(K)

    K <- K / K.transform
    P.temp <- diag(ncol(K)) - matrix(1/ncol(K),ncol(K),ncol(K))
    K.sd.off.diagonal[trait] <- sum(diag(P.temp %*% K %*% P.temp %*% K)) / ncol(K) #sqrt( (ncol(K)-1)*var(K[upper.tri(K)]) + 2*ncol(K)*(ncol(K)-1)*var(K[upper.tri(K)])  )

    K2 <- gwas.obj$kinship2[available.genotypes,available.genotypes]
    K2 <- K2 / KinshipTransform(K2)

    Z.geno              <- createZmatrixForFactor(factor((pheno.list[[trait]])$genotype,levels=available.genotypes))

    pheno.frame <- pheno.list[[trait]]
    pheno.frame <- data.frame(pheno.frame,nested.block.factor=block.factor)
    pheno.frame <- cbind(pheno.frame,matrix(0,nrow(pheno.frame),1))
    names(pheno.frame)[-c(1:4)] <- paste("sim",as.character(trait),sep="")

    if (add.pca) {pheno.frame <- cbind(pheno.frame,gwas.obj$pca[pheno.frame$genotype,])}
    n.pca <- ncol(gwas.obj$pca)

    #######################
    #old

    #if (n.qtl <= gwas.obj$nchr) {
    #  chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
    #} else {
    # chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl,replace=T))
    #}


    #qtl.loc <- rep(0,n.qtl)
    #for (j in 1:n.qtl) {
    #  candidates.j <- candidates[[chr.loc[j]]]
    #  qtl.loc[j] <- sample(candidates.j,size=1)
    #}
    #qtl.loc <- sort(qtl.loc)
    ####################
    #new

    too.much.ld <- T

    while (too.much.ld) {

      if (n.qtl <= gwas.obj$nchr) {
        chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
      } else {
       chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl,replace=T))
      }
      chr.table <- table(chr.loc)

      qtl.loc <- integer(0)
      for (chr in as.integer(names(chr.table))) {qtl.loc <- c(qtl.loc,sample(candidates[[chr]],size=as.integer(chr.table[which(chr == as.integer(names(chr.table)))])))}
      qtl.loc <- sort(qtl.loc)

      qtl.matrix <- t(as.matrix(gwas.obj$markers[qtl.loc,(pheno.list[[trait]])$genotype]))
      #b          <- lm(pheno.frame[,5] ~ qtl.matrix)

      qtl.freqs    <- apply(gwas.obj$markers[qtl.loc,(pheno.list[[trait]])$genotype],1,mean)
      effect.sizes <- sqrt(rel.effects * (sigmaG2 * (n.available-1)/n.available)   /    (qtl.freqs * (1 - qtl.freqs))) #  * gwas.obj$n
      effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

      #explained.phenotypic.variance <- summary(b)$r.squared
      #cat(qtl.matrix); cat(effect.sizes)
      qtl.variance.exact            <- as.numeric(matrix(effect.sizes,nrow=1) %*% cov(qtl.matrix) %*% matrix(effect.sizes,nrow=n.qtl))
      qtl.variance.independent      <- sum(effect.sizes^2 * qtl.freqs*(1-qtl.freqs))

      # cat(max(abs(cor(qtl.matrix)[upper.tri(cor(qtl.matrix))])),'\n')
      if (sum(abs(effect.sizes)) > 0) {
        if (min(qtl.variance.independent,qtl.variance.exact) / max(qtl.variance.independent,qtl.variance.exact) > 0.97) {too.much.ld <- F}
      } else {
        too.much.ld <- F
      }
    }

    ##################

    # check
    # match(colnames(Z.geno),GWAS.obj$plant.names)
    # all(colnames(Z.geno)==colnames(K))
    # all(colnames(Z.geno)==unique(pheno.list[[trait]]$genotype))

    qtl.freqs    <- apply(gwas.obj$markers[qtl.loc,(pheno.list[[trait]])$genotype],1,mean)
    effect.sizes <- sqrt(rel.effects * (sigmaG2 * (n.available-1)/n.available)   /    (qtl.freqs * (1 - qtl.freqs))) #  * gwas.obj$n
    effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

    pos.frame[,trait]    <- qtl.loc
    chr.frame[,trait]    <- chr.loc
    effect.frame[,trait] <- effect.sizes

    g         <- MatrixRoot(K) %*% matrix(rnorm(n.available,sd=sqrt(sigmaG2)))
    g.repl    <- as.numeric(Z.geno %*% g)
    qtl.part  <- as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,colnames(K)]))
    qtl.part.repl  <- as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,(pheno.list[[trait]])$genotype]))
    new.trait <- g.repl + rnorm(n.rep * n.available, sd=sqrt(sigmaE2)) +  qtl.part.repl + as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))

    breeding.values[[trait]] <- data.frame(genotype=colnames(K),value= g + qtl.part,predicted1=rep(NA,ncol(K)),predicted2=rep(NA,ncol(K)),row.names=colnames(K))

    pheno.frame[,5]         <- new.trait
    # moved to below:
    #    pheno.list.out[[trait]] <- pheno.frame

    qtl.matrix <- t(as.matrix(gwas.obj$markers[qtl.loc,(pheno.list[[trait]])$genotype]))
    b          <- lm(pheno.frame[,5] ~ qtl.matrix)

    explained.phenotypic.variance[trait] <- summary(b)$r.squared
    qtl.variance.exact[trait]            <- as.numeric(matrix(effect.sizes,nrow=1) %*% cov(qtl.matrix) %*% matrix(effect.sizes,nrow=n.qtl))
    qtl.variance.independent[trait]      <- sum(effect.sizes^2 * qtl.freqs*(1-qtl.freqs))

    #######################################################
    # setting some observations to missing

    if (n.corrupt.geno > 0) {
      candidate.genotypes <- names(which(table(pheno.frame$genotype[!is.na(pheno.frame[,5])]) > n.left)) # n.rep - n.corrupt.blocks
      stopifnot(length(candidate.genotypes) >= n.corrupt.geno)
      corrupt.genos <- sample(candidate.genotypes,size=n.corrupt.geno) #corrupt.genos <- sample(GWAS.obj$plant.names,size=n.corrupt.geno)
      for (cg in corrupt.genos) {
        candidate.plants <- which(pheno.frame$genotype==cg & !is.na(pheno.frame[,5]))
        remaining.ones <- sample(candidate.plants,size=sample(1:n.left,size=1)) # remaining.ones <- sample(candidates,size=1)
        pheno.frame[setdiff(candidate.plants,remaining.ones),5] <- NA
      }
    }
    pheno.frame <- pheno.frame[!is.na(pheno.frame[,5]),]
    pheno.list.out[[trait]] <- pheno.frame

    temp.pheno <- aggregate(pheno.frame[,5], by=list(pheno.frame$genotype),FUN=mean,na.rm =T)
    names(temp.pheno) <- c('genotype','sim')
    rownames(temp.pheno) <- temp.pheno$genotype
    if (add.pca) {temp.pheno <- cbind(temp.pheno,gwas.obj$pca[temp.pheno$genotype,])}


    ######################## data-simulation is now complete ; estimation starts from here

    anova.h2.estimate[trait]       <- Heritability(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)[1]
    anova.h2.estimate2[trait]      <- Heritability3(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)[1]

    conf.int.left.broad.sense[trait] <- Heritability3(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)$conf.int.left
    conf.int.right.broad.sense[trait] <- Heritability3(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)$conf.int.right

    #########

    # load(file='test_gs.RData')

    #n.rep.vector.temp <- table(pheno.frame$genotype)[pheno.frame$genotype]
    n.rep.vector.temp <- as.integer(table(pheno.frame$genotype))
    #average.number.of.replicates.temp <- 1 / mean(1/n.rep.vector.temp)
    average.number.of.replicates.temp <- ( sum(n.rep.vector.temp) - sum(n.rep.vector.temp^2) / sum(n.rep.vector.temp) ) / ( length(n.rep.vector.temp) - 1 )

    if (use.asreml) {

      if (line.heritability) {av.temp <- average.number.of.replicates.temp} else {av.temp <- 1}

      if (add.pca) {

        reml.obj1B <- HeritabilityMixedModel(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype,
                                             covariates=pheno.frame[,-(1:5)],K=K2[unique(pheno.frame$genotype),unique(pheno.frame$genotype)],
                                             no.BLUPs=no.BLUPs,line.heritability=line.heritability,average.number.of.replicates=av.temp)
        rm(fixed.formula,random.formula)

      } else {
        reml.obj1B <- HeritabilityMixedModel(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype,
                                             K=K2[unique(pheno.frame$genotype),unique(pheno.frame$genotype)],
                                             no.BLUPs=no.BLUPs,line.heritability=line.heritability,average.number.of.replicates=av.temp)
        rm(fixed.formula,random.formula)
      }

  	  #reml.obj1B<- HeritabilityMixedModel(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype,K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)])
      conf.int1.left.1stage[trait] <- reml.obj1B$conf.int1[1]
      conf.int1.right.1stage[trait] <- reml.obj1B$conf.int1[2]
      conf.int2.left.1stage[trait] <- reml.obj1B$conf.int2[1]
      conf.int2.right.1stage[trait] <- reml.obj1B$conf.int2[2]

      mixed.model.h2.estimate.1stageB[trait] <- reml.obj1B$h2#reml.obj1B$vg / (reml.obj1B$vg + reml.obj1B$ve)

    } else {

      if (add.pca) {
    	  #reml.obj1 <- emma.REMLE(as.numeric(pheno.frame[,5]),as.matrix(cbind(as.matrix(rep(1,n.available*n.rep)),pheno.frame[,-(1:5)])),K2[pheno.frame$genotype,pheno.frame$genotype])
        reml.obj1 <- emma.REMLE(as.numeric(pheno.frame[,5]),as.matrix(cbind(as.matrix(rep(1,nrow(pheno.frame))),pheno.frame[,-(1:5)])),K2[pheno.frame$genotype,pheno.frame$genotype])
      } else {
    	  #reml.obj1 <- emma.REMLE(pheno.frame[,5],as.matrix(rep(1,n.available*n.rep)),K2[pheno.frame$genotype,pheno.frame$genotype])
    	  reml.obj1 <- emma.REMLE(pheno.frame[,5],as.matrix(rep(1,nrow(pheno.frame))),K2[pheno.frame$genotype,pheno.frame$genotype])
      }
      #n.rep.vector <- as.integer(table(pheno.frame$genotype))
      #n.rep.vector <- rep(n.rep.vector,times=n.rep.vector)
      #reml.obj1B <- emma.REMLE(pheno.frame[,5],as.matrix(n.rep.vector),K[pheno.frame$genotype,pheno.frame$genotype])

      if (line.heritability) {
        mixed.model.h2.estimate.1stage[trait]  <- reml.obj1$vg / (reml.obj1$vg + reml.obj1$ve / average.number.of.replicates.temp) # mean(table(pheno.frame$genotype))
      } else {
        mixed.model.h2.estimate.1stage[trait]  <- reml.obj1$vg / (reml.obj1$vg + reml.obj1$ve )
      }
    }

    if (!no.BLUPs) {

      if (use.asreml) {

        breeding.values.MSE1[trait]          <- mean((reml.obj1B$blup.frame[as.character(breeding.values[[trait]]$genotype),2] - breeding.values[[trait]]$value)^2)
        #@#new
        breeding.values[[trait]]$predicted1  <- reml.obj1B$blup.frame[as.character(breeding.values[[trait]]$genotype),2]

      } else {

        temp.blup.frame <- data.frame(genotype=unique(pheno.frame$genotype),predicted.value=rep(NA,length(unique(pheno.frame$genotype))))
        row.names(temp.blup.frame) <- temp.blup.frame$genotype
        temp.blup.frame$genotype <- as.character(temp.blup.frame$genotype)

        delta.temp <- reml.obj1$vg / reml.obj1$ve # var.comp[1] / var.comp[2]
        K.temp     <- K2#[pheno.frame$genotype,pheno.frame$genotype]
        V          <- delta.temp * K.temp[pheno.frame$genotype,pheno.frame$genotype] + diag(nrow(pheno.frame))
        W          <- solve(V)
        mu         <- sum(W %*% as.matrix(pheno.frame[,5])) / sum(W)

        Z.temp <- createZmatrixForFactor(f.vector=as.character(pheno.frame$genotype),factor.name="",ordered.factor=F)
        Z.temp <- Z.temp[,colnames(K.temp)]
        temp.blup.frame[colnames(Z.temp),]$predicted.value <- delta.temp * K.temp %*% t(Z.temp) %*% W %*% as.matrix(pheno.frame[,5]-mu) + mu
        breeding.values.MSE1[trait] <- mean((temp.blup.frame[as.character(breeding.values[[trait]]$genotype),2] - breeding.values[[trait]]$value)^2)
        #@#new
        breeding.values[[trait]]$predicted1  <- temp.blup.frame[as.character(breeding.values[[trait]]$genotype),2]

      }
    }

    if (use.asreml) {

      if (line.heritability) {av.temp <- 1} else {av.temp <- average.number.of.replicates.temp}

      if (add.pca) {
        reml.obj2B <- HeritabilityMixedModel(data.vector=temp.pheno$sim,geno.vector=temp.pheno$genotype,covariates=temp.pheno[,-(1:2)],
                                             K=K2[temp.pheno$genotype,temp.pheno$genotype],no.BLUPs=no.BLUPs,line.heritability=line.heritability,
                                             average.number.of.replicates=av.temp)
        rm(fixed.formula,random.formula)
      } else {
     	  reml.obj2B <- HeritabilityMixedModel(data.vector=temp.pheno$sim,geno.vector=temp.pheno$genotype,
                                             K=K2[temp.pheno$genotype,temp.pheno$genotype],no.BLUPs=no.BLUPs,line.heritability=line.heritability,
                                             average.number.of.replicates=av.temp)
        if (MSE.method) {
       	  reml.obj3B <- HeritabilityMixedModelWithFixedSigmaE2(data.vector=temp.pheno$sim,geno.vector=temp.pheno$genotype,
                                               K=K2[temp.pheno$genotype,temp.pheno$genotype],no.BLUPs=no.BLUPs,line.heritability=line.heritability,
                                               average.number.of.replicates=av.temp,sigmaE2=Heritability3(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype)[3])
        }

        rm(fixed.formula,random.formula)
      }

  	  #reml.obj1B<- HeritabilityMixedModel(data.vector=pheno.frame[,5],geno.vector=pheno.frame$genotype,K=K[unique(pheno.frame$genotype),unique(pheno.frame$genotype)])
      conf.int1.left.2stage[trait] <- reml.obj2B$conf.int1[1]
      conf.int1.right.2stage[trait] <- reml.obj2B$conf.int1[2]

      conf.int2.left.2stage[trait] <- reml.obj2B$conf.int2[1]
      conf.int2.right.2stage[trait] <- reml.obj2B$conf.int2[2]

      mixed.model.h2.estimate.2stageB[trait] <- reml.obj2B$h2#reml.obj2$vg / (reml.obj2$vg + reml.obj2$ve)
      #
      if (MSE.method) {
        mixed.model.h2.estimate.2stageBextra[trait] <- reml.obj3B$h2
        mixed.model.sg.estimate.2stageBextra[trait] <- reml.obj3B$vg
        mixed.model.se.estimate.2stageBextra[trait] <- reml.obj3B$ve
      }
    } else {
      if (add.pca) {
     	  reml.obj2 <- emma.REMLE(temp.pheno$sim,as.matrix(cbind(as.matrix(rep(1,n.available)),temp.pheno[,-(1:2)])),K2[temp.pheno$genotype,temp.pheno$genotype])
      } else {
     	  reml.obj2 <- emma.REMLE(temp.pheno$sim,as.matrix(rep(1,n.available)),K2[temp.pheno$genotype,temp.pheno$genotype])
      }
      if (line.heritability) {
        mixed.model.h2.estimate.2stage[trait]  <- reml.obj2$vg / (reml.obj2$vg + reml.obj2$ve)
      } else {
        mixed.model.h2.estimate.2stage[trait]  <- reml.obj2$vg / (reml.obj2$vg + reml.obj2$ve * average.number.of.replicates.temp)
      }
    }

    if (!no.BLUPs) {
      if (use.asreml) {
        breeding.values.MSE2[trait] <- mean((reml.obj2B$blup.frame[as.character(breeding.values[[trait]]$genotype),2] - breeding.values[[trait]]$value)^2)
        #@#new
        breeding.values[[trait]]$predicted2  <- reml.obj2B$blup.frame[as.character(breeding.values[[trait]]$genotype),2]

      } else {

        temp.blup.frame <- data.frame(genotype=unique(pheno.frame$genotype),predicted.value=rep(NA,length(unique(pheno.frame$genotype))))
        row.names(temp.blup.frame) <- temp.blup.frame$genotype
        temp.blup.frame$genotype <- as.character(temp.blup.frame$genotype)

        delta.temp <- reml.obj2$vg / reml.obj2$ve # var.comp[1] / var.comp[2]
        K.temp <- K2#[temp.frame$genotype,temp.frame$genotype]
        V <- delta.temp * K.temp[temp.pheno$genotype,temp.pheno$genotype] + diag(nrow(temp.pheno))
        W <- solve(V)
        mu <- sum(W %*% as.matrix(temp.pheno[,2])) / sum(W)

        Z.temp <- createZmatrixForFactor(f.vector=as.character(temp.pheno$genotype),factor.name="",ordered.factor=F)
        Z.temp <- Z.temp[,colnames(K.temp)]
        temp.blup.frame[colnames(Z.temp),]$predicted.value <- delta.temp * K.temp %*% t(Z.temp) %*% W %*% as.matrix(temp.pheno[,2]-mu) + mu
        breeding.values.MSE2[trait] <- mean((temp.blup.frame[as.character(breeding.values[[trait]]$genotype),2] - breeding.values[[trait]]$value)^2)
        #@#new
        breeding.values[[trait]]$predicted2  <- temp.blup.frame[as.character(breeding.values[[trait]]$genotype),2]

      }
    }
    #mixed.model.h2.estimate.1stage[trait] <- reml.obj1$vg / (reml.obj1$vg + reml.obj1$ve)
    #mixed.model.h2.estimate.2stage[trait] <- reml.obj2$vg / (reml.obj2$vg + reml.obj2$ve)
    #mixed.model.h2.estimate.1stage[trait] <- reml.obj1$vg / (reml.obj1$vg + reml.obj1$ve / n.rep)

    if (use.asreml) {
      sigmaG2.estimate.1[trait] <- reml.obj1B$vg
      sigmaG2.estimate.2[trait] <- reml.obj2B$vg
      sigmaE2.estimate.1[trait] <- reml.obj1B$ve
      sigmaE2.estimate.2[trait] <- reml.obj2B$ve
    } else {
      sigmaG2.estimate.1[trait] <- reml.obj1$vg
      sigmaG2.estimate.2[trait] <- reml.obj2$vg
      sigmaE2.estimate.1[trait] <- reml.obj1$ve
      sigmaE2.estimate.2[trait] <- reml.obj2$ve
    }
    #narrow.h2.estimate.1stage[trait]  <- reml.obj1$vg / var(pheno.frame[,5])
    #narrow.h2.estimate.2stage[trait]  <- reml.obj2$vg / var(temp.pheno$sim)


    if (n.pred > 0) {

      GS[[trait]] <- data.frame(genotype=sample(setdiff(gwas.obj$plant.names,colnames(K)),size=n.pred),value=rep(NA,n.pred),
                                predicted1=rep(NA,n.pred),predicted2=rep(NA,n.pred))
      GS[[trait]]$genotype <- as.character(GS[[trait]]$genotype)
      row.names(GS[[trait]]) <- GS[[trait]]$genotype

      K.pred <- gwas.obj$kinship[GS[[trait]]$genotype,GS[[trait]]$genotype]
      K.pred <- K.pred / K.transform

      K.pred.obs <- gwas.obj$kinship[GS[[trait]]$genotype,colnames(K)]
      K.pred.obs <- K.pred.obs / K.transform

      #GS[[trait]]$value <- as.numeric(K.pred.obs %*% solve(K) %*% g)  + as.numeric(MatrixRoot(K.pred.obs %*% solve(K) %*% t(K.pred.obs)) %*% matrix(rnorm(n.pred,sd=sqrt(sigmaG2))))
      #GS[[trait]]$value <- as.numeric(K.pred.obs %*% solve(K) %*% g)  + as.numeric(t(chol((K.pred.obs %*% solve(K) %*% t(K.pred.obs)))) %*% matrix(rnorm(n.pred,sd=sqrt(sigmaG2))))
      GS[[trait]]$value <- as.numeric(K.pred.obs %*% solve(K) %*% g)  + as.numeric(MatrixRoot(K.pred - K.pred.obs %*% solve(K) %*% t(K.pred.obs)) %*% matrix(rnorm(n.pred,sd=sqrt(sigmaG2))))
      GS[[trait]]$value <- GS[[trait]]$value  + as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,GS[[trait]]$genotype]))

      W.1 <- solve(sigmaG2.estimate.1[trait] * K[pheno.frame$genotype,pheno.frame$genotype] + sigmaE2.estimate.1[trait] * diag(nrow(pheno.frame)))
      W.2 <- solve(sigmaG2.estimate.2[trait] * K[temp.pheno$genotype,temp.pheno$genotype] + sigmaE2.estimate.2[trait] * diag(nrow(temp.pheno)))

      mu.hat.1 <- sum(W.1 %*% as.matrix(pheno.frame[,5])) / sum(W.1)
      mu.hat.2 <- sum(W.2 %*% as.matrix(temp.pheno[,2])) / sum(W.2)

      GS[[trait]]$predicted1 <- mu.hat.1 + sigmaG2.estimate.1[trait] * as.numeric(K.pred.obs[,pheno.frame$genotype] %*% W.1 %*% as.matrix(pheno.frame[,5] - mu.hat.1))
      GS[[trait]]$predicted2 <- mu.hat.2 + sigmaG2.estimate.2[trait] * as.numeric(K.pred.obs[,temp.pheno$genotype] %*% W.2 %*% as.matrix(temp.pheno[,2] - mu.hat.2))
      #GS[[trait]]$predicted2 <- sigmaG2.estimate.2[trait] * as.numeric(K.pred.obs %*% solve(sigmaG2.estimate.2[trait] * K + sigmaE2.estimate.2[trait] * diag(ncol(K))) %*% as.matrix(temp.pheno[colnames(K),2]))

      GS.MSE1[trait] <- mean((GS[[trait]]$value-GS[[trait]]$predicted1)^2)
      GS.MSE2[trait] <- mean((GS[[trait]]$value-GS[[trait]]$predicted2)^2)

      check.predictions <- F
      #save.image(file='test_GS.RData')
      #  reml.pred$predictions$vcov

      if (check.predictions) {

        ############
        # 2-stage : check with asreml
        #library(asreml)
        qw <- data.frame(genotype=c(temp.pheno$genotype,GS[[trait]]$genotype),value=c(temp.pheno[,2],rep(NA,length(GS[[trait]]$genotype))))
        row.names(qw) <- qw$genotype

        fixed.formula <<- as.formula('value ~ 1')
        random.formula <<- as.formula("~ giv(genotype,var=T)")
        qw$genotype <- as.character(qw$genotype)
        asreml.temp <<- MakeKinshipAsreml(GWAS.obj$kinship[qw$genotype,qw$genotype])
        reml.obj <- asreml(fixed.formula,data=qw,random = random.formula,ginverse = list(genotype=asreml.temp))

        reml.pred <- predict(reml.obj, classify='genotype',data=qw)
        blup.frame <- reml.pred$predictions$pvals
        blup.frame$genotype <- as.character(blup.frame$genotype)
        row.names(blup.frame) <- blup.frame$genotype

        print(GS[[trait]][,]$predicted2 - blup.frame[GS[[trait]]$genotype,]$predicted.value)
        #blup.frame[order(blup.frame$standard.error),][1:210,]
        #setdiff(blup.frame[order(blup.frame$standard.error),][201:250,1],GS[[1]][,1])

        average.pred.error2[trait] <- mean(reml.pred$predictions$pvals$standard.error,na.rm=T)

        ##################

        N <- nrow(reml.obj2B$blup.frame)
        sp <- sigmaE2.estimate.2[trait] * solve((sigmaE2.estimate.2[trait] / sigmaG2.estimate.2[trait]) * solve(K[temp.pheno$genotype,temp.pheno$genotype]) + diag(N) - matrix(1/N,N,N))
        diag(sp)
        km <- rbind(matrix(1,1,N),diag(N))
        Cmatrix     <- matrix(0,N+1,N+1)
        Cmatrix[1,] <- cbind(N,matrix(1,1,N))
        Cmatrix[-1,1]  <- matrix(1,N,1)
        Cmatrix[-1,-1] <- (sigmaE2.estimate.2[trait] / sigmaG2.estimate.2[trait]) * solve(K[temp.pheno$genotype,temp.pheno$genotype]) + diag(N)# - matrix(1/N,N,N))
        CmatrixInv <- solve(Cmatrix)
        sp2 <- sigmaE2.estimate.2[trait] * t(km) %*% CmatrixInv %*% km
        diag(sp2)

        cor(reml.obj2B$blup.frame$standard.error,diag(sp))
        summary( reml.obj2B$blup.frame$standard.error - diag(sp))
        all(reml.obj2B$blup.frame$genotype==colnames(sp))


        sigmaG2.estimate.2[trait] * diag(K.pred - K.pred.obs %*% solve(K) %*% t(K.pred.obs))
        diag(K.pred.obs %*% solve(K) %*% sp  %*% t(K.pred.obs %*% solve(K)))

        ##############
        # 1-stage : check with asreml
        qw <- data.frame(genotype=c(pheno.frame$genotype,GS[[trait]]$genotype),value=c(pheno.frame[,5],rep(NA,length(GS[[trait]]$genotype))))
        #row.names(qw) <- qw$genotype

        fixed.formula <<- as.formula('value ~ 1')
        random.formula <<- as.formula("~ giv(genotype,var=T)")
        qw$genotype <- as.character(qw$genotype)
        asreml.temp <<- MakeKinshipAsreml(GWAS.obj$kinship[unique(qw$genotype),unique(qw$genotype)])
        reml.obj <- asreml(fixed.formula,data=qw,random = random.formula,ginverse = list(genotype=asreml.temp))

        reml.pred <- predict(reml.obj, classify='genotype',data=qw)
        blup.frame <- reml.pred$predictions$pvals
        blup.frame$genotype <- as.character(blup.frame$genotype)
        row.names(blup.frame) <- blup.frame$genotype

        print(GS[[trait]][,]$predicted1 - blup.frame[GS[[trait]]$genotype,]$predicted.value)
        cor(GS[[trait]][,]$predicted1 , blup.frame[GS[[trait]]$genotype,]$predicted.value )

        average.pred.error1[trait] <- mean(reml.pred$predictions$pvals$standard.error,na.rm=T)
      }

    }

  }

  sigma.G.total <- sigmaG2 * (1 + sum(rel.effects) )
  #h2 <- sigma.G.total / (sigma.G.total + sigmaE2 * (n.available-1) / (n.available) )

  simulated.h2 <- sigma.G.total / (sigma.G.total + sigmaE2 * (n.available-1) / (n.available * n.rep ) )
  simulated.sigma.G.total <- sigma.G.total
  simulated.sigmaG2 <- sigmaG2
  simulated.QTL.variance <- sigmaG2 * sum(rel.effects)

return(list(pheno.list=pheno.list.out,pos.frame=pos.frame,chr.frame=chr.frame,effect.frame=effect.frame,simulated.h2=simulated.h2,
            simulated.sigma.G.total=simulated.sigma.G.total,simulated.sigmaG2=simulated.sigmaG2,simulated.QTL.variance=simulated.QTL.variance,
            rel.effects=rel.effects,sigmaG2=sigmaG2,sigmaE2=sigmaE2,sigmaB2=sigmaB2,
            explained.phenotypic.variance=explained.phenotypic.variance,qtl.variance.exact=qtl.variance.exact,qtl.variance.independent=qtl.variance.independent,
            h2.1stage=mixed.model.h2.estimate.1stage,h2.1stageB=mixed.model.h2.estimate.1stageB,
            h2.2stage=mixed.model.h2.estimate.2stage,h2.2stageB=mixed.model.h2.estimate.2stageB,
            h2.broad.sense=unlist(anova.h2.estimate),h2.broad.sense2=unlist(anova.h2.estimate2),
            h2.estimate.2stageBextra=mixed.model.h2.estimate.2stageBextra,
            sg.estimate.2stageBextra=mixed.model.sg.estimate.2stageBextra,
            se.estimate.2stageBextra=mixed.model.se.estimate.2stageBextra,
            sigmaG2.estimate.1=sigmaG2.estimate.1,sigmaG2.estimate.2=sigmaG2.estimate.2,sigmaE2.estimate.1=sigmaE2.estimate.1,sigmaE2.estimate.2=sigmaE2.estimate.2,
            conf.int1.left.1stage=conf.int1.left.1stage,conf.int1.right.1stage=conf.int1.right.1stage,
            conf.int1.left.2stage=conf.int1.left.2stage,conf.int1.right.2stage=conf.int1.right.2stage,
            conf.int2.left.1stage=conf.int2.left.1stage,conf.int2.right.1stage=conf.int2.right.1stage,
            conf.int2.left.2stage=conf.int2.left.2stage,conf.int2.right.2stage=conf.int2.right.2stage,
            conf.int3.left.1stage=conf.int3.left.1stage,conf.int3.right.1stage=conf.int3.right.1stage,
            conf.int3.left.2stage=conf.int3.left.2stage,conf.int3.right.2stage=conf.int3.right.2stage,
            conf.int.left.broad.sense=unlist(conf.int.left.broad.sense), conf.int.right.broad.sense=unlist(conf.int.right.broad.sense),
            K.sd.off.diagonal=K.sd.off.diagonal,
            breeding.values=breeding.values,breeding.values.MSE1=breeding.values.MSE1,breeding.values.MSE2=breeding.values.MSE2,
            GS=GS,GS.MSE1=GS.MSE1,GS.MSE2=GS.MSE2,
            average.pred.error1=average.pred.error1,average.pred.error2=average.pred.error2)) #h2.anova2=unlist(anova.h2.estimate2),
            #narrow.h2.estimate.1stage=narrow.h2.estimate.1stage,narrow.h2.estimate.2stage=narrow.h2.estimate.2stage
}






SimulatePhenoGivenGenoAlhpaDesign_with_column_effect <- function(gwas.obj,n.sim=10,rel.effects=c(0.01,0.02,0.05,0.1),
                                                                 sigmaA2=150,sigmaE2=65,sigmaB2=200,sigmaC2=100,different.names=T) {
# This is a variation on the function SimulatePhenoGivenGenoAlhpaDesign : here, we add an extra column effect
#
# new version: qtl-locations are drawn again for each simulation
# simulation of a phenotype given real genotypic GWAS data.
# An (already generated) alpha design is assumed
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
#   columns 2 and 3 of gwas.obj$pheno must be replicate and block
# * n.sim : number of simulated data sets
# * different.names : when TRUE, the simulated traits will have col-names sim1, sim2,...
#                    otherwise they will all have col-name pheno
# OUTPUT :
# *
# *
# gwas.obj=GWAS.obj;n.sim=10;rel.effects=c(0.001,0.003);sigmaA2=0.1;sigmaE2=1;sigmaB2=1;sigmaC2=1;different.names=T
  #for the time being, assume that length(rel.effects) is at most the number of chromosomes

  n.qtl <- length(rel.effects)
  n.rep <- unique(table(gwas.obj$pheno$genotype))# N.B. assumes equal numbers of replicates !!!                [1] #table(gwas.obj$genotype)
  available.genotypes <- unique(gwas.obj$pheno$genotype)
  n.available         <- length(available.genotypes)

  mafs <- list()
  candidates <- list()
  for (j in gwas.obj$chromosomes) {  # N.B. GWAS.obj$chromosomes should be labeled 1,2,3...
    mafs[[j]] <- apply(gwas.obj$markers[gwas.obj$map$chromosome==j,],1,mean)
    candidates[[j]] <- ((1:gwas.obj$N)[gwas.obj$map$chromosome==j])[which(mafs[[j]] > .1 & mafs[[j]] < .9)]
  }

  pos.frame <- as.data.frame(matrix(0,n.qtl,n.sim))
  names(pos.frame) <- paste("sim",as.character(1:n.sim),sep="")
  chr.frame <- pos.frame
  effect.frame <- pos.frame

  #z <- c(rep(c(rep(1,n.rep),rep(0,n.rep*gwas.obj$n)),gwas.obj$n-1),rep(1,n.rep))
  #Z.geno <- matrix(z,ncol=gwas.obj$n)
  #Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=gwas.obj$plant.names))
  Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=available.genotypes))

  Z.block <- createZmatrixForIncompleteBlockDesign(rep.vector=gwas.obj$pheno[,2],block.vector=gwas.obj$pheno[,3])
  block.factor  <- Z.block$nested.block.factor
  Z.block <- Z.block$Zmatrix

  Z.column <- createZmatrixForFactor(factor(gwas.obj$pheno$column,levels=levels(gwas.obj$pheno$column)))

  pheno.frame        <- gwas.obj$pheno
  pheno.frame <- data.frame(pheno.frame,nested.block.factor=block.factor)
  pheno.frame <- cbind(pheno.frame,matrix(0,nrow(pheno.frame),n.sim))

  if (different.names) {
    names(pheno.frame)[-c(1:5)] <- paste("sim",as.character(1:n.sim),sep="")
  } else {
    names(pheno.frame)[-c(1:5)] <- rep("pheno",n.sim)
  }

  for (trait in 1:n.sim) {
    chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
    qtl.loc <- rep(0,n.qtl)
    for (j in 1:n.qtl) {
      #mafs.j <- mafs[[chr.loc[j]]]
      candidates.j <- candidates[[chr.loc[j]]]
      qtl.loc[j] <- sample(candidates.j,size=1)
    }
    qtl.freqs <- apply(gwas.obj$markers[qtl.loc,],1,mean)
    effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship[available.genotypes,available.genotypes]) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    #effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

    pos.frame[,trait] <- qtl.loc
    chr.frame[,trait] <- chr.loc
    effect.frame[,trait] <- effect.sizes

    # na.omit(unique(gwas.obj$pheno$genotype)[match(gwas.obj$plant.names,unique(gwas.obj$pheno$genotype))])
    g <- MatrixRoot(gwas.obj$kinship[available.genotypes,available.genotypes]) %*% matrix(rnorm(n.available,sd=sqrt(sigmaA2)))
    g.repl <- as.numeric(Z.geno %*% g)
    new.trait <- g.repl + rnorm(n.rep * n.available, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype])) + as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2)))) + as.numeric(Z.column %*% matrix(rnorm(ncol(Z.column),sd=sqrt(sigmaC2))))
    pheno.frame[,5+trait] <- new.trait

    # old code
    #g <- MatrixRoot(gwas.obj$kinship) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))
    #g.repl <- as.numeric(Z.geno %*% g)
    #new.trait <- g.repl + rnorm(n.rep * gwas.obj$n, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype]))+ as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))
    #pheno.frame[,4+trait] <- new.trait
    ##             #+  rep(as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,])),each=n.rep)

  }

return(list(pheno.frame=pheno.frame,pos.frame=pos.frame,chr.frame=chr.frame,effect.frame=effect.frame))
}




SimulatePhenoGivenGenoAlhpaDesign <- function(gwas.obj,n.sim=10,rel.effects=c(0.01,0.02,0.05,0.1),sigmaA2=150,sigmaE2=65,sigmaB2=200,different.names=T) {
# new version: qtl-locations are drawn again for each simulation
# simulation of a phenotype given real genotypic GWAS data.
# An (already generated) alpha design is assumed
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
#   columns 2 and 3 of gwas.obj$pheno must be replicate and block
# * n.sim : number of simulated data sets
# * different.names : when TRUE, the simulated traits will have col-names sim1, sim2,...
#                    otherwise they will all have col-name pheno
# OUTPUT :
# *
# *
# gwas.obj=GWAS.obj;n.sim=10;rel.effects=c(0.001,0.003);sigmaA2=0.1;sigmaE2=1;sigmaB2=200;different.names=T
  #for the time being, assume that length(rel.effects) is at most the number of chromosomes

  n.qtl <- length(rel.effects)
  n.rep <- unique(table(gwas.obj$pheno$genotype))# N.B. assumes equal numbers of replicates !!!                [1] #table(gwas.obj$genotype)
  available.genotypes <- unique(gwas.obj$pheno$genotype)
  n.available         <- length(available.genotypes)

  mafs <- list()
  candidates <- list()
  for (j in gwas.obj$chromosomes) {  # N.B. GWAS.obj$chromosomes should be labeled 1,2,3...
    mafs[[j]] <- apply(gwas.obj$markers[gwas.obj$map$chromosome==j,],1,mean)
    candidates[[j]] <- ((1:gwas.obj$N)[gwas.obj$map$chromosome==j])[which(mafs[[j]] > .1 & mafs[[j]] < .9)]
  }

  pos.frame <- as.data.frame(matrix(0,n.qtl,n.sim))
  names(pos.frame) <- paste("sim",as.character(1:n.sim),sep="")
  chr.frame <- pos.frame
  effect.frame <- pos.frame

  #z <- c(rep(c(rep(1,n.rep),rep(0,n.rep*gwas.obj$n)),gwas.obj$n-1),rep(1,n.rep))
  #Z.geno <- matrix(z,ncol=gwas.obj$n)
  #Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=gwas.obj$plant.names))
  Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=available.genotypes))


  Z.block <- createZmatrixForIncompleteBlockDesign(rep.vector=gwas.obj$pheno[,2],block.vector=gwas.obj$pheno[,3])
  block.factor  <- Z.block$nested.block.factor
  Z.block <- Z.block$Zmatrix

  pheno.frame        <- gwas.obj$pheno
  pheno.frame <- data.frame(pheno.frame,nested.block.factor=block.factor)
  pheno.frame <- cbind(pheno.frame,matrix(0,nrow(pheno.frame),n.sim))

  if (different.names) {
    names(pheno.frame)[-c(1:4)] <- paste("sim",as.character(1:n.sim),sep="")
  } else {
    names(pheno.frame)[-c(1:4)] <- rep("pheno",n.sim)
  }

  for (trait in 1:n.sim) {
    chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
    qtl.loc <- rep(0,n.qtl)
    for (j in 1:n.qtl) {
      #mafs.j <- mafs[[chr.loc[j]]]
      candidates.j <- candidates[[chr.loc[j]]]
      qtl.loc[j] <- sample(candidates.j,size=1)
    }
    qtl.freqs <- apply(gwas.obj$markers[qtl.loc,],1,mean)
    effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship[available.genotypes,available.genotypes]) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    #effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

    pos.frame[,trait] <- qtl.loc
    chr.frame[,trait] <- chr.loc
    effect.frame[,trait] <- effect.sizes

    # na.omit(unique(gwas.obj$pheno$genotype)[match(gwas.obj$plant.names,unique(gwas.obj$pheno$genotype))])
    g <- MatrixRoot(gwas.obj$kinship[available.genotypes,available.genotypes]) %*% matrix(rnorm(n.available,sd=sqrt(sigmaA2)))
    g.repl <- as.numeric(Z.geno %*% g)
    new.trait <- g.repl + rnorm(n.rep * n.available, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype]))+ as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))
    pheno.frame[,4+trait] <- new.trait

    # old code
    #g <- MatrixRoot(gwas.obj$kinship) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))
    #g.repl <- as.numeric(Z.geno %*% g)
    #new.trait <- g.repl + rnorm(n.rep * gwas.obj$n, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype]))+ as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))
    #pheno.frame[,4+trait] <- new.trait
    ##             #+  rep(as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,])),each=n.rep)

  }

return(list(pheno.frame=pheno.frame,pos.frame=pos.frame,chr.frame=chr.frame,effect.frame=effect.frame))
}

SimulatePhenoGivenGenoAlhpaDesign_with_column_effect <- function(gwas.obj,n.sim=10,rel.effects=c(0.01,0.02,0.05,0.1),
                                                                 sigmaA2=150,sigmaE2=65,sigmaB2=200,sigmaC2=100,different.names=T) {
# This is a variation on the function SimulatePhenoGivenGenoAlhpaDesign : here, we add an extra column effect
#
# new version: qtl-locations are drawn again for each simulation
# simulation of a phenotype given real genotypic GWAS data.
# An (already generated) alpha design is assumed
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
#   columns 2 and 3 of gwas.obj$pheno must be replicate and block
# * n.sim : number of simulated data sets
# * different.names : when TRUE, the simulated traits will have col-names sim1, sim2,...
#                    otherwise they will all have col-name pheno
# OUTPUT :
# *
# *
# gwas.obj=GWAS.obj;n.sim=10;rel.effects=c(0.001,0.003);sigmaA2=0.1;sigmaE2=1;sigmaB2=1;sigmaC2=1;different.names=T
  #for the time being, assume that length(rel.effects) is at most the number of chromosomes

  n.qtl <- length(rel.effects)
  n.rep <- unique(table(gwas.obj$pheno$genotype))# N.B. assumes equal numbers of replicates !!!                [1] #table(gwas.obj$genotype)
  available.genotypes <- unique(gwas.obj$pheno$genotype)
  n.available         <- length(available.genotypes)

  mafs <- list()
  candidates <- list()
  for (j in gwas.obj$chromosomes) {  # N.B. GWAS.obj$chromosomes should be labeled 1,2,3...
    mafs[[j]] <- apply(gwas.obj$markers[gwas.obj$map$chromosome==j,],1,mean)
    candidates[[j]] <- ((1:gwas.obj$N)[gwas.obj$map$chromosome==j])[which(mafs[[j]] > .1 & mafs[[j]] < .9)]
  }

  pos.frame <- as.data.frame(matrix(0,n.qtl,n.sim))
  names(pos.frame) <- paste("sim",as.character(1:n.sim),sep="")
  chr.frame <- pos.frame
  effect.frame <- pos.frame

  #z <- c(rep(c(rep(1,n.rep),rep(0,n.rep*gwas.obj$n)),gwas.obj$n-1),rep(1,n.rep))
  #Z.geno <- matrix(z,ncol=gwas.obj$n)
  #Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=gwas.obj$plant.names))
  Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=available.genotypes))

  Z.block <- createZmatrixForIncompleteBlockDesign(rep.vector=gwas.obj$pheno[,2],block.vector=gwas.obj$pheno[,3])
  block.factor  <- Z.block$nested.block.factor
  Z.block <- Z.block$Zmatrix

  Z.column <- createZmatrixForFactor(factor(gwas.obj$pheno$column,levels=levels(gwas.obj$pheno$column)))

  pheno.frame        <- gwas.obj$pheno
  pheno.frame <- data.frame(pheno.frame,nested.block.factor=block.factor)
  pheno.frame <- cbind(pheno.frame,matrix(0,nrow(pheno.frame),n.sim))

  if (different.names) {
    names(pheno.frame)[-c(1:5)] <- paste("sim",as.character(1:n.sim),sep="")
  } else {
    names(pheno.frame)[-c(1:5)] <- rep("pheno",n.sim)
  }

  for (trait in 1:n.sim) {
    chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
    qtl.loc <- rep(0,n.qtl)
    for (j in 1:n.qtl) {
      #mafs.j <- mafs[[chr.loc[j]]]
      candidates.j <- candidates[[chr.loc[j]]]
      qtl.loc[j] <- sample(candidates.j,size=1)
    }
    qtl.freqs <- apply(gwas.obj$markers[qtl.loc,],1,mean)
    effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship[available.genotypes,available.genotypes]) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    #effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

    pos.frame[,trait] <- qtl.loc
    chr.frame[,trait] <- chr.loc
    effect.frame[,trait] <- effect.sizes

    # na.omit(unique(gwas.obj$pheno$genotype)[match(gwas.obj$plant.names,unique(gwas.obj$pheno$genotype))])
    g <- MatrixRoot(gwas.obj$kinship[available.genotypes,available.genotypes]) %*% matrix(rnorm(n.available,sd=sqrt(sigmaA2)))
    g.repl <- as.numeric(Z.geno %*% g)
    new.trait <- g.repl + rnorm(n.rep * n.available, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype])) + as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2)))) + as.numeric(Z.column %*% matrix(rnorm(ncol(Z.column),sd=sqrt(sigmaC2))))
    pheno.frame[,5+trait] <- new.trait

    # old code
    #g <- MatrixRoot(gwas.obj$kinship) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))
    #g.repl <- as.numeric(Z.geno %*% g)
    #new.trait <- g.repl + rnorm(n.rep * gwas.obj$n, sd=sqrt(sigmaE2)) +  as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,gwas.obj$pheno$genotype]))+ as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))
    #pheno.frame[,4+trait] <- new.trait
    ##             #+  rep(as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,])),each=n.rep)

  }

return(list(pheno.frame=pheno.frame,pos.frame=pos.frame,chr.frame=chr.frame,effect.frame=effect.frame))
}


SimulatePhenoGivenTesterDesign <- function(gwas.obj,n.sim=10,rel.effects=c(0.01,0.02,0.05,0.1),sigmaA2=150,sigmaE2=65,sigmaB2=200) { # ,rep.variance=0
# new version: qtl-locations are drawn again for each simulation
# simulation of a phenotype given real genotypic GWAS data.
# An (already generated) alpha design is assumed
# INPUT :
# * gwas.obj : a GWAS.obj, whose markers field is a p x n.ind  matrix or data.frame of 0/1 marker-scores
#   columns 2 and 3 of gwas.obj$pheno must be replicate and block
# * n.sim : number of simulated data sets
# * different.names : when TRUE, the simulated traits will have col-names sim1, sim2,...
#                    otherwise they will all have col-name pheno
# OUTPUT :
# *
# *

#gwas.obj =GWAS.obj

  n.qtl <- length(rel.effects)
  #n.rep <- unique(table(gwas.obj$pheno$genotype))# N.B. assumes equal numbers of replicates !!!                [1] #table(gwas.obj$genotype)
  n.rep <- 1
  mafs <- list()
  candidates <- list()
  for (j in gwas.obj$chromosomes) {  # N.B. GWAS.obj$chromosomes should be labeled 1,2,3...
    mafs[[j]] <- apply(gwas.obj$markers[gwas.obj$map$chromosome==j,],1,mean)
    candidates[[j]] <- ((1:gwas.obj$N)[gwas.obj$map$chromosome==j])[which(mafs[[j]] > .1 & mafs[[j]] < .9)]
  }


  pos.frame <- as.data.frame(matrix(0,n.qtl,n.sim))
  names(pos.frame) <- paste("sim",as.character(1:n.sim),sep="")
  chr.frame <- pos.frame
  effect.frame <- pos.frame
  Z.geno <- createZmatrixForFactor(factor(gwas.obj$pheno$genotype,levels=gwas.obj$plant.names))

  Z.block <- createZmatrixForIncompleteBlockDesign(rep.vector=gwas.obj$pheno[,2],block.vector=gwas.obj$pheno[,3])
  block.factor  <- Z.block$nested.block.factor
  Z.block <- Z.block$Zmatrix

  pheno.frame        <- gwas.obj$pheno
  pheno.frame <- data.frame(pheno.frame,nested.block.factor=block.factor)
  pheno.frame <- cbind(pheno.frame,matrix(0,nrow(pheno.frame),n.sim))


  #names(pheno.frame)[-c(1:4)] <- rep("pheno",n.sim)
  names(pheno.frame)[-c(1:4)] <- paste("sim",as.character(1:n.sim),sep="")

  for (trait in 1:n.sim) {
    chr.loc <- sort(sample(gwas.obj$chromosomes,size=n.qtl))
    qtl.loc <- rep(0,n.qtl)
    for (j in 1:n.qtl) {
      #mafs.j <- mafs[[chr.loc[j]]]
      candidates.j <- candidates[[chr.loc[j]]]
      qtl.loc[j] <- sample(candidates.j,size=1)
    }
    qtl.freqs <- apply(gwas.obj$markers[qtl.loc,],1,mean)
    effect.sizes <- sqrt(rel.effects * sigmaA2 * KinshipTransformUnscaled(gwas.obj$kinship) / (qtl.freqs * (1 - qtl.freqs) * gwas.obj$n))
    effect.sizes <- effect.sizes * sample(c(-1,1),replace=T,size=n.qtl)

    pos.frame[,trait] <- qtl.loc
    chr.frame[,trait] <- chr.loc
    effect.frame[,trait] <- effect.sizes

    #g <- MatrixRoot(gwas.obj$kinship) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))
    g <- t(chol(gwas.obj$kinship)) %*% matrix(rnorm(gwas.obj$n,sd=sqrt(sigmaA2)))

    g <- g - mean(g)
    g.rep <- as.numeric(Z.geno %*% g)

    ########### !!!!!!!!!!!!!!
    # added in the following line : as.character
    qtl.signal.rep <- as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,as.character(gwas.obj$pheno$genotype)]))
    #qtl.signal <- as.numeric(gwas.obj$markers[qtl.loc,])

    new.trait <- g.rep + rnorm(nrow(Z.geno), sd=sqrt(sigmaE2)) + qtl.signal.rep + as.numeric(Z.block %*% matrix(rnorm(ncol(Z.block),sd=sqrt(sigmaB2))))
    pheno.frame[,4+trait] <- new.trait
                 #+  rep(as.numeric( matrix(effect.sizes,nrow=1) %*% as.matrix(gwas.obj$markers[qtl.loc,])),each=n.rep)

  }

return(list(pheno.frame=pheno.frame,pos.frame=pos.frame,chr.frame=chr.frame,effect.frame=effect.frame,Gsim=g.rep,qtl.signal.rep=qtl.signal.rep))
}


plotPhenoInAlphaDesign <- function(y,rep.vector,block.vector,x.rep,y.rep,x.block,y.block,x.plot,y.plot) {
# OBJECTIVE : construct a matrix from a phenotypic vector y. The elements of y are positioned according to the field layout; the output of this
#             function can be used as input in an image command.
#             An alpha design is assumed
# where the (x,y) coordinates are determined by
#
# INPUT :
# * y : vector of phenotypic values
# * rep.vector : vector with the replicates the observations y were taken
# * block.vector : vector with the blocks the observations y were taken
#   Remarks : 1) the blocks are assumed to be numbered 1....B WITHIN every replicate
#             2) y, rep.vector and block.vector are typically columns in a gwas-object (in GWAS.obj$pheno)
# * x.rep, y.rep : the number of complete replicates in the horizontal (x) and vertical (y) direction.
#                    In total there are x.rep*y.rep complete replicates
# * x.block, y.block : the number of blocks within each complete replicate; in the horizontal (x) and vertical (y) direction.
#                    In total there are x.block*y.block  blocks within each complete replicate
# * x.plot, y.plot : the number of plots within each block in the horizontal (x) and vertical (y) direction.
#                    In total there are x.plot*y.plot  plots within each block
#  The length of  y needs to be x.plot*x.block*x.rep*y.plot*y.block*y.rep
#
# OUTPUT : a matrix, which can directly be used to create a 2-D field map : using image(function_output)
#
# EXAMPLE :
# qw<- plotPhenoInAlphaDesign(y=GWAS.obj$pheno[,6],rep.vector=GWAS.obj$pheno$rep,block.vector=GWAS.obj$pheno$block,x.rep=4,y.rep=2,x.block=2,y.block=7,x.plot=5,y.plot=5)
# image(qw)

  rep.vector <- as.integer(rep.vector)
  block.vector <- as.integer(block.vector)
  n.x <- x.rep*x.block*x.plot
  n.y <- y.rep*y.block*y.plot
  result <- matrix(0,n.y,n.x)
  xx <- (floor((rep.vector-1) %% x.rep) * x.block + floor((block.vector-1) %% x.block) ) * x.plot #block.vector
  yy <- (floor((rep.vector-1)/x.rep) * y.block + floor((block.vector-1)/x.block) ) * y.plot
  xx <- xx + rep(rep(1:x.plot,y.plot),n.x*n.y/(x.plot*y.plot))
  yy <- yy + rep(rep(1:y.plot,each=x.plot),n.x*n.y/(x.plot*y.plot))
  result[cbind(yy,xx)] <- y
  result <- as.matrix(result)
return(t(result))
}

WriteKinshipFile <- function(M,plant.names,file.name) {
# OBJECTIVE :
# write a kinship matrix to a file
# INPUT :
# plant.names : a vector of genotype names (as character)
# M : the kinship matrix
  cat(paste(plant.names,collapse=","),"\n",file=file.name)
  write.table(M,file=file.name,append=T,quote=T,row.names=F,col.names=F,sep=",")
}
WriteKinshipFile <- function(M,plant.names,file.name) {
# OBJECTIVE :
# write a kinship matrix to a file
# INPUT :
# plant.names : a vector of genotype names (as character)
# M : the kinship matrix
  M <- as.data.frame(M)
  names(M) <- as.character(plant.names)
  write.table(M,file=file.name,quote=T,row.names=F,col.names=T,sep=",")
}
# WriteKinshipFile(diag(GWAS.obj$n),GWAS.obj$plant.names,"LFN_identity.csv")


###############################################################################################################

CreateGroupsBasedOnGwasOutput <- function(pvalues,gwas.obj,window.size,max.number=length(pvalues)) {
# OBJECTIVE : given a (physical) map and a vector of p-values, construct a vector of markers (more precisely, marker numbers) that are most significant, and that are not in LD.
#   First, the most significant marker is selected, and a region around this marker is excluded.
#   Subsequent markers are found similarly, by selecting the most significant markers outside the 'LD-regions' of the markers that have already been selected.
#   The process stops if either these LD-regions completely cover the complete genome, or if there no markers left with p-value smaller than a threshold (default min.p=0.1)
# INPUT :
#   gwas.obj : a GWAS-object, with gwas.obj$N markers
#   pvalues : a vector of pvalues fo length gwas.obj$N, obtained for gwas.obj
#   window.size : this number indicates how many base -pairs on both sides of a selected marker are excluded
# OUTPUT :
#   the marker numbers, i.e. element numbers in the vector pvalues
#
#browser()
  group.begin.end.matrix <- data.frame(start.snp=integer(0),end.snp=integer(0)) # matrix(0,0,2)
  stop.search <- F
  available.markers <- rep(T,gwas.obj$N)
  pvalues[is.na(pvalues)] <- 1
  selected.markers <- integer(0)
  m <- 0
  while (!stop.search) {
    smallest.p <- min(pvalues[available.markers])
    best.marker <- which.min(pvalues[available.markers])[1]  # N.B if several markers have the same (minimal) pvalue, we choose the first one
    best.marker <- ((1:gwas.obj$N)[available.markers])[best.marker]
    selected.markers <- c(selected.markers,best.marker)
    latest.region <- sort(intersect(GetSNPsInRegion(gwas.obj,best.marker,region.size=window.size),(1:gwas.obj$N)[available.markers]))
    available.markers[c(latest.region,best.marker)] <- F
    new.block.start <- min(c(latest.region,best.marker))
    new.block.end   <- max(c(latest.region,best.marker))
    group.begin.end.matrix <- rbind(group.begin.end.matrix,c(new.block.start,new.block.end))
    #cat(latest.region,"\n\n")
    #cat(c(new.block.start,new.block.end),"\n\n\n\n")
    if (sum(available.markers)==0) {stop.search <- T}
    m <- m + 1
    if (m >= max.number) {stop.search <- T}
  }
  names(group.begin.end.matrix) <- c("start.snp","end.snp")
  group.begin.end.matrix <- group.begin.end.matrix[order(group.begin.end.matrix$start.snp),]
return(group.begin.end.matrix)
}
#group.matrix <- CreateGroupsBasedOnGwasOutput(GWA.result$pvalue,GWAS.obj,400000)

AddGroupsToGwasObj <- function(group.matrix,gwas.obj) {
  K          <- nrow(group.matrix)
  if (length(unique(group.matrix[,2]-group.matrix[,1]))==1) {
    group.list <- data.frame(apply(group.matrix,1,function(x){(x[1]):(x[2])}))
  } else {
    group.list <- apply(group.matrix,1,function(x){(x[1]):(x[2])})
  }
  group.vector <- rep(0,gwas.obj$N)
  for (i in 1:K) {group.vector[group.matrix[i,1]:group.matrix[i,2]] <- i}
  names(group.list) <- 1:K
  gwas.obj$groups <- list(group.list=group.list,n.group=K,group.vector=group.vector,position.matrix=group.matrix)
return(gwas.obj)
}
#GWAS.obj <- AddGroupsToGwasObj(group.matrix,gwas.obj=GWAS.obj)

REMLlikelihoodDelta <- function(gwas.obj,tr.n,h2,REML=T) {
    # For a given vector of heritability (h2) values, this function
    # calculates the REML log-likelihood for the 'empty - model', i.e. without markers,
    # for the trait in column tr.n of gwas.obj$pheno
    #
    # gwas.obj=GWAS.obj;tr.n=31;h2=0.8;REML=T

    gwas.obj$pheno$genotype <- as.character(gwas.obj$pheno$genotype)
    kinship.reduced  <- MakeScanGlsKinship(1,0,gwas.obj$kinship,plant.names=gwas.obj$plant.names,gwas.obj$pheno,tr.n=tr.n)
    Knr <- kinship.reduced / KinshipTransform(kinship.reduced)
    #Knr <- kinship.reduced
      #
      # The following lines follow the derivation on p. 3 in the supplement of Lippert et a;. 2012 Nat. Genet.
      # N.B. This is for a model with only an intercept and no snp, so d=1
      #
      evk <- eigen(Knr)
      U <- evk[[2]]
      v <- evk[[1]]
      xbar <- apply(U,2,sum)
      ybar <- as.numeric(t(U) %*% as.matrix(gwas.obj$pheno[!is.na(gwas.obj$pheno[,tr.n]),tr.n]))


   if (length(h2)==1) {
      delta <- (1-h2) / h2
      #delta <- h2 / (1-h2)
      sbar <- 1/(v+delta)
      n <- length(sbar)
      betahat <- sum(sbar * xbar * ybar) / sum(sbar * xbar * xbar)
      sigmaG2 <- sum(sbar * (ybar - betahat * xbar)^2) / (n-1)  # supplement of Lippert et a;. 2012 Nat. Genet. : p. 10 , second display
      if (REML) {
        LLdelta <- -0.5 * ((n-1)*log(2*pi) + n -log(n) - sum(log(sbar)) + (n-1) * log(sigmaG2)  + log(sum(U %*% diag(sbar) %*% t(U))))
      } else {
        LLdelta <- -0.5 * (n*log(2*pi) + n - sum(log(sbar)) + n * log(sigmaG2) )
      }
    } else {
      LLdelta <- rep(NA,length(h2))
      for (i in 1:length(h2)) {
        delta <- (1-h2[i]) / h2[i]
        sbar <- 1/(v+delta)
        n <- length(sbar)
        betahat <- sum(sbar * xbar * ybar) / sum(sbar * xbar * xbar)
        sigmaG2 <- sum(sbar * (ybar - betahat * xbar)^2) / (n-1)  # supplement of Lippert et a;. 2012 Nat. Genet. : p. 10 , second display
        if (REML) {
          LLdelta[i] <- -0.5 * ((n-1)*log(2*pi) + n -log(n) - sum(log(sbar)) + (n-1) * log(sigmaG2)  + log(sum(U %*% diag(sbar) %*% t(U))))
        } else {
          LLdelta[i] <- -0.5 * (n*log(2*pi) + n - sum(log(sbar)) + n * log(sigmaG2) )
        }
      }
    }
return(LLdelta)
}

PlotREMLlikelihoodDelta <- function(gwas.obj,tr.n,main.title=names(gwas.obj$pheno)[tr.n],n.points=100,file.name=NULL) {
  h2 <- (1:n.points)/(n.points+1)
  if (is.null(file.name)) {
    plot(h2,REMLlikelihoodDelta(gwas.obj=gwas.obj,tr.n=tr.n,h2=h2),type="l",ylab="log-likelihood",main=main.title)
  } else {
    jpeg(quality=100,file=file.name)
    plot(h2,REMLlikelihoodDelta(gwas.obj=gwas.obj,tr.n=tr.n,h2=h2),type="l",ylab="log-likelihood",main=main.title)
    dev.off()
  }
}


fastGLS <-function(Y,X,Sigma) {  # gwas.obj,
#test: Y=Y.temp;X=X.temp;cofs=Z
  # check for missing values, and class of Y
  n1     <- length(Y)
  M1     <- solve(chol(Sigma))
  Y_t1   <- crossprod(M1,Y)             # pre-multiply the phenotype (Y1_) with t(M1)
  int_t1 <- crossprod(M1,(rep(1,n1)))   # pre-multiply the intercept with t(M1)

  if (ncol(X)==1) { # for extra robustness, distinguish....
    X_t1 <- crossprod(M1,matrix(as.numeric(X)))           # pre-multiply the snp-matrix with t(M1)
  } else {
    X_t1 <- crossprod(M1,X)           # pre-multiply the snp-matrix with t(M1)
  }

  RSS_env<-rep(sum(lsfit(int_t1,Y_t1,intercept = FALSE)$residuals^2),ncol(X))
  R1_full<-apply(X_t1,2,function(x){sum(lsfit(cbind(int_t1,x),Y_t1,intercept = FALSE)$residuals^2)})

  F_1<-((RSS_env-R1_full)/1)/(R1_full/(n1-2))
  pval_Y1<-pf(F_1,1,(n1-2),lower.tail=FALSE)
  rm(X_t1,Y_t1,M1)
  #out_models1<-data.frame(SNP=names(pval_Y1),Pval_Y1=pval_Y1)
#return(out_models1)
return(pval_Y1)
}

fastGLS_with_effect_sizes <- function(Y,X,Sigma) {
# Y=Y.temp;X=X.temp
  # check for missing values, and class of Y
  n1     <- length(Y)
  M1     <- solve(chol(Sigma))
  Y_t1   <- crossprod(M1,Y)             # pre-multiply the phenotype (Y1_) with t(M1)
  int_t1 <- crossprod(M1,(rep(1,n1)))   # pre-multiply the intercept with t(M1)

  if (ncol(X)==1) { # for extra robustness, distinguish....
    X_t1 <- crossprod(M1,matrix(as.numeric(X)))           # pre-multiply the snp-matrix with t(M1)
  } else {
    X_t1 <- crossprod(M1,X)           # pre-multiply the snp-matrix with t(M1)
  }

  ###########    matrix cookbook, 3.2.6 Rank-1 update of inverse of inner product
  a  <- 1 / (as.numeric(t(int_t1) %*% int_t1))
  VV <- X_t1 * X_t1
  vv <- apply(VV,2,sum)
  vX <- as.numeric(t(int_t1) %*% X_t1)
  nn <- 1 / (vv - a * vX^2)
  XtXinv2ndRows <- cbind(-1* a * vX * nn,nn)
  Xty <- cbind(rep(as.numeric(t(int_t1) %*% Y_t1),length(nn)),as.numeric(t(Y_t1) %*% X_t1))
  beta_vector <- XtXinv2ndRows[,1] * Xty[,1] + XtXinv2ndRows[,2] * Xty[,2] #apply(XtXinv2ndRows * Xty,2,sum)
  ##########
  # check, for the first snp
  #G <- matrix(0,nrow(X.temp),2)
  #G[,1] <- as.numeric(int_t1)
  #G[,2] <- as.numeric(X_t1[,2])
  #solve(t(G) %*% G) %*% t(G) %*% Y_t1
  ###############
  #solve(t(G) %*% G)
  #1 / (vv[1] / a - (vX[1])^2*a)
  #- a * vX[1] / (vv[1] - (vX[1])^2*a)
  ##################

  RSS_env <- rep(sum(lsfit(int_t1,Y_t1,intercept = FALSE)$residuals^2),ncol(X))
  R1_full <- apply(X_t1,2,function(x){sum(lsfit(cbind(int_t1,x),Y_t1,intercept = FALSE)$residuals^2)})

  F_1     <- ((RSS_env-R1_full)/1)/(R1_full/(n1-2))
  pval_Y1 <- pf(F_1,1,(n1-2),lower.tail=FALSE)
  rm(X_t1,Y_t1,M1)

  # the R_LR^2 statistic from G. Sun et al 2010, heredity
  R_LR_2  <- 1 - exp(-(1 / n1)*(RSS_env-R1_full))

  #out_models1<-data.frame(SNP=names(pval_Y1),Pval_Y1=pval_Y1)
#return(out_models1)
return(list(pvalue=pval_Y1,beta=beta_vector,beta_se=sqrt(XtXinv2ndRows[,2]),R_LR_2 = R_LR_2))
}

fastGLS_cof_with_effect_sizes <-function(Y,X,Sigma,cofs,nbchunks=10) {
# Y=Y.temp;X=X.temp;cofs=Z; nbchunks=10#Sigma=Sigma
# code taken from mlmm_cof.r (Segura et al (2012) Nature Gen.) , and simplified
#
# OBJECTIVE : compute p-values for the GLS F-test as in emma-x, Fast-LMM or scan_GLS
#
#
# INPUT :
# Y : vector of phenotypic values (type numeric or integer, with names)
# X : n x p matrix of marker-scores  (n being the number of individuals, p the number of markers). Of type matrix
# Sigma : the covariance matrix
# cofs  : the n x d matrix  of covariates (NOT including the intercept)
# nbchunks : the number of parts (chunks) to in which the calculations are split
#
# Important:
#
# - no missing values allowed in any of the input
# - everything should be in the same order, and the rownames of Y,X and Sigma should be the same (as well as the column-names of Sigma, and the names of Y)
# - the colnames of X should be the marker names (X is of class matrix, so use colnames and not col.names)
# - there may be replicates of the same genotype, but these should then be all contained in X
#
# OUTPUT :
# * a vector of pvalues
#
# Questions/to check:
# *
#
#test: Y=Y.temp;X=X.temp;cofs=createZmatrixForFactor(GWAS.obj$pheno[!is.na(GWAS.obj$pheno[,tr.n]),rep.col]);nbchunks<- 10

# Y=Y.temp;X=X.temp;cofs=Z;Sigma=Sigma;nbchunks=10
  df1      <- 1

  n        <- length(Y)
  m        <- ncol(X)

  Xo       <- rep(1,n)
  fix_cofs <- cbind(Xo,cofs)
  n.cov.extra <- ncol(cofs)


  COF      <- fix_cofs

  M        <- solve(chol(Sigma))
  Y_t      <- crossprod(M,Y)

  cof_fwd_t <- crossprod(M,COF)
  fwd_lm   <- summary(lm(Y_t~0+cof_fwd_t))
  Res_H0   <- fwd_lm$residuals
  Q_       <- qr.Q(qr(cof_fwd_t))

  ###############

  RSS<-list()

  for (j in 1:(nbchunks-1)) {
  #j=1
    X_t      <- crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    RSS[[j]] <- apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
    rm(X_t)
  }
  X_t <- crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),X[,-(1:((j)*round(m/nbchunks)))])
  RSS[[nbchunks]] <- apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})

  rm(X_t,j)

  ################

  n1     <- length(Y)
  M1     <- solve(chol(Sigma))
  Y_t1   <- crossprod(M1,Y)             # pre-multiply the phenotype (Y1_) with t(M1)
  #int_t1 <- crossprod(M1,(rep(1,n1)))   # pre-multiply the intercept with t(M1)

  COF_t      <- crossprod(M1,COF)

  if (ncol(X)==1) { # for extra robustness, distinguish....
    X_t1 <- crossprod(M1,matrix(as.numeric(X)))           # pre-multiply the snp-matrix with t(M1)
  } else {
    X_t1 <- crossprod(M1,X)           # pre-multiply the snp-matrix with t(M1)
  }

  A  <- solve(t(COF_t) %*% COF_t)

  VV <- X_t1 * X_t1
  vv <- apply(VV,2,sum)

  vX <- t(COF_t) %*% X_t1     # c x N

  nn <- 1 / (vv - apply(vX * (A %*% vX),2,sum))

  XtXinvLastRows <- cbind(- nn * t(vX) %*% A,nn)

  Xty <- cbind(matrix(rep((as.numeric(t(COF_t) %*% Y_t1)),length(nn)),byrow=T,ncol=n.cov.extra+1),as.numeric(t(Y_t1) %*% X_t1))

  beta_vector <- apply(XtXinvLastRows[,1:(1+n.cov.extra)] * Xty[,1:(1+n.cov.extra)],1,sum) + XtXinvLastRows[,2+n.cov.extra] * Xty[,2+n.cov.extra]

  ##############

  RSSf    <- unlist(RSS)
  RSS_H0  <- sum(Res_H0^2)
  df2     <- n - df1 - ncol(fix_cofs)#- ncol(cof_fwd[[1]])
  Ftest   <- (rep(RSS_H0,length(RSSf))/RSSf-1)*df2/df1
  pval    <- pf(Ftest,df1,df2,lower.tail=FALSE)

  # the R_LR^2 statistic from G. Sun et al 2010, heredity
  R_LR_2  <- 1 - exp(-(1 / n1)*(rep(RSS_H0,length(RSSf))-RSSf))

return(list(pvalue=pval,beta=beta_vector,beta_se=sqrt(nn),R_LR_2 = R_LR_2))
}



fastGLS_cof <-function(Y,X,Sigma,cofs,nbchunks=10) {
# code taken from mlmm_cof.r (Segura et al (2012) Nature Gen.) , and simplified
#
# OBJECTIVE : compute p-values for the GLS F-test as in emma-x, Fast-LMM or scan_GLS
#
#
# INPUT :
# Y : vector of phenotypic values (type numeric or integer, with names)
# X : n x p matrix of marker-scores  (n being the number of individuals, p the number of markers). Of type matrix
# Sigma : the covariance matrix
# cofs  : the n x d matrix  of covariates (NOT including the intercept)
# nbchunks : the number of parts (chunks) to in which the calculations are split
#
# Important:
#
# - no missing values allowed in any of the input
# - everything should be in the same order, and the rownames of Y,X and Sigma should be the same (as well as the column-names of Sigma, and the names of Y)
# - the colnames of X should be the marker names (X is of class matrix, so use colnames and not col.names)
# - there may be replicates of the same genotype, but these should then be all contained in X
#
# OUTPUT :
# * a vector of pvalues
#
# Questions/to check:
# *
#
#test: Y=Y.temp;X=X.temp;cofs=createZmatrixForFactor(GWAS.obj$pheno[!is.na(GWAS.obj$pheno[,tr.n]),rep.col]);nbchunks<- 10
  df1<-1

  n<-length(Y)
  m<-ncol(X)

  Xo<-rep(1,n)
  fix_cofs<-cbind(Xo,cofs)
  #COF <- cbind(fix_cofs,cof_fwd[[1]])
  COF <- fix_cofs

  M<-solve(chol(Sigma))
  Y_t<-crossprod(M,Y)

  #cof_fwd_t <- crossprod(M,cbind(fix_cofs,cof_fwd[[1]]))
  cof_fwd_t <- crossprod(M,COF)
  fwd_lm <- summary(lm(Y_t~0+cof_fwd_t))
  Res_H0 <- fwd_lm$residuals
  Q_ <- qr.Q(qr(cof_fwd_t))

  RSS<-list()
  for (j in 1:(nbchunks-1)) {
    X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
    rm(X_t)
  }
  X_t<-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),X[,-(1:((j)*round(m/nbchunks)))])
  RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})

  rm(X_t,j)

  RSSf    <- unlist(RSS)
  RSS_H0  <- sum(Res_H0^2)
  df2     <- n - df1 - ncol(fix_cofs)#- ncol(cof_fwd[[1]])
  Ftest   <- (rep(RSS_H0,length(RSSf))/RSSf-1)*df2/df1
  pval    <- pf(Ftest,df1,df2,lower.tail=FALSE)

return(pval)
}

fastGLS_cof_single_marker <-function(Y,X,Sigma,cofs) {
# code taken from mlmm_cof.r (Segura et al (2012) Nature Gen.) , and simplified
# Compared to fastGLS_cof : here we do it for a single marker
#
# OBJECTIVE : compute p-values for the GLS F-test as in emma-x, Fast-LMM or scan_GLS
#
#
# INPUT :
# Y : vector of phenotypic values (type numeric or integer, with names)
# X : n x 1 matrix of marker-scores  (n being the number of individuals, p the number of markers). Of type matrix
# Sigma : the covariance matrix
# cofs  : the n x d matrix  of covariates (NOT including the intercept)
#
# Important:
#
# - no missing values allowed in any of the input
# - everything should be in the same order, and the rownames of Y,X and Sigma should be the same (as well as the column-names of Sigma, and the names of Y)
# - the colnames of X should be the marker names (X is of class matrix, so use colnames and not col.names)
# - there may be replicates of the same genotype, but these should then be all contained in X
#
# OUTPUT :
# * a vector of pvalues
#
# Questions/to check:
# *
#
#test: Y=Y.temp;X=X.temp;cofs=createZmatrixForFactor(GWAS.obj$pheno[!is.na(GWAS.obj$pheno[,tr.n]),rep.col])

  df1<-1

  n<-length(Y)
  m<-ncol(X)

  Xo<-rep(1,n)
  fix_cofs<-cbind(Xo,cofs)
  #COF <- cbind(fix_cofs,cof_fwd[[1]])
  COF <- fix_cofs

  M <- solve(chol(Sigma))
  Y_t <- crossprod(M,Y)

  #cof_fwd_t <- crossprod(M,cbind(fix_cofs,cof_fwd[[1]]))
  cof_fwd_t <- crossprod(M,COF)
  fwd_lm <- summary(lm(Y_t~0+cof_fwd_t))
  Res_H0 <- fwd_lm$residuals
  Q_ <- qr.Q(qr(cof_fwd_t))

  X_t <-crossprod(M %*% (diag(n)-tcrossprod(Q_,Q_)),matrix(as.numeric(X)))

  RSS <-apply(X_t,2,function(x){sum(lsfit(x,Res_H0,intercept = FALSE)$residuals^2)})
  rm(X_t)

  RSSf    <- unlist(RSS)
  RSS_H0  <- sum(Res_H0^2)
  df2     <- n - df1 - ncol(fix_cofs)#- ncol(cof_fwd[[1]])
  Ftest   <- (rep(RSS_H0,length(RSSf))/RSSf-1)*df2/df1
  pval    <- pf(Ftest,df1,df2,lower.tail=FALSE)

return(pval)
}


AddPca <- function(pheno.frame,vars,prefix="",scaled.version=T) {
# N.B. only use on genotypic means
# missing values are removed

# OBJECTIVE : add the principle components scores to a phenotypic data.frame
# *
# INPUT :
# * pheno.frame : a data-frame of phenotypic data, with one column per trait
# * vars : vector of column numbers corresponding to pheno.frame; indicates over which variables the PCA is to be performed
# * prefix : the new PC-variables have names paste(prefix,"pca",as.character(1:ncol(d$scores)),sep="")
# * scaled.version : PCA on standardized variables ?
# OUTPUT :
# * the same phenotypic data.frame, but with the principal components added

  a <-pheno.frame[,vars]
  a <- a[apply(a,1,function(x){sum(is.na(x))})==0,]
  #b <- princomp(x=a) # PCA without scaling
  d <- princomp(x=a,cor=scaled.version)
  new.traits <- d$scores
  pheno.frame <- cbind(pheno.frame,matrix(NA,nrow(pheno.frame),ncol(d$scores)))
  pheno.frame[row.names(new.traits),(ncol(pheno.frame)-ncol(d$scores)+1):(ncol(pheno.frame))] <- new.traits
  names(pheno.frame)[(ncol(pheno.frame)-ncol(d$scores)+1):(ncol(pheno.frame))] <- paste(prefix,"pca",as.character(1:ncol(d$scores)),sep="")
return(pheno.frame)
}
# OBJECTIVE :
# *
# *
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *


GwasObjToMtmmObj <- function(gwas.obj) {

  mtmm.obj <- list()
  mtmm.obj$n <- gwas.obj$n
  mtmm.obj$N <- gwas.obj$N
  mtmm.obj$Y <- gwas.obj$pheno
  mtmm.obj$K <- gwas.obj$kinship

  mtmm.obj$X    <- matrix(rep(as.integer(0),gwas.obj$N*gwas.obj$n),ncol=gwas.obj$N)
  blocks <- DefineBlocks(indices=1:gwas.obj$N,block.size=20000)

  for (b in 1:length(blocks)) {
    mtmm.obj$X[,blocks[[b]]] <- t(gwas.obj$markers[blocks[[b]],])
    gc()
  }

  colnames(mtmm.obj$X) <- paste(as.character(gwas.obj$map$chromosome),"- ",as.character(gwas.obj$map$position),sep="")
  rownames(mtmm.obj$X) <- gwas.obj$plant.names

return(mtmm.obj)
}

ProcessCycdesignOutput <- function(file.name,n.block.x=1,rep.layout="ver",fill.reps.by.col=F,turn.upside.down=F) {
# INPUT :
# - file.name : file name of the cyc-design output (type .txt ; no .des or .html !)
# - n.block.x : number of blocks next to each other
# - rep.layout :  options : "hor", "ver". With hor(izontal), the replicates are assumed to be next to each other
# - fill.reps.by.col : determines whether (within each replicate) the blocks will be numbered column by column, or row by row. See the example below.
# - turn.upside.down : after making the layout, should the field be put 'upside down' ? (put to TRUE if you want to start in the bottom left corner)
#
# OUTPUT : a list with the following elements
# - n.rep : the number of replicates (determined from the information in the file file.name)
# - n.geno : the number of genotypes (determined from the information in the file file.name)
# - rep.list : list of length n.rep ; each element is a data-frame, indicating which genotype-numbers are in which block.
#              (number of columns = the number of incomplete blocks within a replicate; number of rows= number of plots within in an incomplete block)
# - layout.list : list of length n.rep ; each element is a data-frame, indicating the spatial layout of the genotype-numbers.
#                 row and column names of the data-frames depend on the parameters  rep.layout and  turn.upside.down
# - rep.layout : the rep.layout value that was given as input
# - turn.upside.down : the turn.upside.down value that was given as input
# - spreadsheet : data-frame , with one row for each plot. Contains empty columns  Order_Nr Code_ID Hybrid_ID

#    example:
#
#    replicate plot block row column Order_Nr Code_ID Hybrid_ID Design_ID
#1           1    1     1   1      1       NA      NA        NA       248
#2           1    2     1   1      2       NA      NA        NA        27
#3           1    3     1   1      3       NA      NA        NA       163
#4           1    4     1   1      4       NA      NA        NA       146
#5           1    5     1   1      5       NA      NA        NA       175
#6           1    6     1   1      6       NA      NA        NA        55
#
#
# COMMENTS :
#
# assumes that blocks are laid out horizontally (next to one another), i.e. the blocks consist of (1 times block.size) plots,
# where block.size is determined from the cyc-design output
# The function automatically detects the number of blocks that should be underneath one another (within each replicate) : n.block.y
# Also the number of genotypes (treatments) is determined from the cyc-design output : n.geno

#Example :
#
#Suppose that the blocks are 1 x 7 plots, and each replicate is laid out as follows:
#
#block1  block2  ... block9
#block10 block11 ... block18
#block19 block20 ... block27
#block28 block29 ... block36
#
#Then fill.reps.by.col=F and n.block.x=9, and n.block.y, n.geno should be respectively 4, 252

#When, for the same design, fill.reps.by.col=T, you get
#block1 block5 ... block33
#block2 block6 ... block34
#block3 block7 ... block35
#block4 block8 ... block36

# turn.upside.down=T, the last layout will change to
#block4 block8 ... block36
#block3 block7 ... block35
#block2 block6 ... block34
#block1 block5 ... block33
# N.B also the order of the rep's will be changed, if rep.layout="ver" !


# file.name='nerac_wd.txt';rep.layout="ver";n.block.x=9;fill.reps.by.col=T
  if (!(rep.layout %in% c('ver','hor')) ) {rep.layout <- 'ver'}

  a <- readChar(file.name, file.info(file.name)$size)
  b <- unlist(strsplit(a,split="plot"))[-1]

  n.rep      <- length(b)
  block.size <- length(unlist(strsplit(b[1],split="\n")))-5


  rep.list <- list()

  for (i in 1:n.rep) {
    write.table(b[i],file="qwerty.txt",col.names=F)
    #d <- read.table('qwerty.txt')
    d <- readLines('qwerty.txt',n=100000)
    write.table(d[2:(block.size+1)],file="asdf.txt",col.names=F,quote=F,row.names=F)
    design <- read.table("asdf.txt")
    design <- design[,-(1:2)]
    names(design) <- paste("block",1:ncol(design),sep="")
    rep.list[[i]] <- design
  }

  n.geno     <-  ncol(rep.list[[1]])*nrow(rep.list[[1]])
  n.block.per.rep <- n.geno / block.size
  n.col <- n.block.x * block.size
  n.block <- n.geno / block.size
  n.block.y <- n.block / n.block.x

  layout.list <- list()
  spreadsheet <- data.frame(replicate=rep(1:n.rep,each=n.geno),plot=rep(1:n.geno,n.rep),block=rep(1:n.block,each=block.size),
                            row=rep(NA,n.rep*n.geno),column=rep(NA,n.rep*n.geno),Order_Nr=rep(NA,n.rep*n.geno),Code_ID=rep(NA,n.rep*n.geno),Hybrid_ID=rep(NA,n.rep*n.geno),
                            #Variety_ID=rep(NA,n.rep*n.geno),Accession_ID=rep(NA,n.rep*n.geno),
                            Design_ID=rep(NA,n.rep*n.geno))

  for (i in 1:n.rep) {
    # old row names :
    #if (n.block.x==1) {
    #  rowNames <- paste('rep',i,'_block',1:n.block.per.rep,sep='')
    #} else {
    #  rowNames <- paste('rep',i,'_block_',(0:(n.block.y-1))*n.block.x+1,"_to_",(0:(n.block.y-1))*n.block.x+n.block.x,sep='')
    #}
    qw <- as.numeric(as.matrix(rep.list[[i]]))

    layout.list[[i]] <- as.data.frame(matrix(qw,ncol=n.col,byrow=T))
    spreadsheet$Design_ID[(i-1)*n.geno + 1:n.geno] <- qw # (n.block.x * block.size)

    if (fill.reps.by.col & n.block.x>1) {
      for (j in 1:n.block.x) {
        layout.list[[i]][,1:block.size + (j-1)* block.size] <- t(as.matrix(rep.list[[i]])[,1:n.block.y + (j-1)*(n.block.y)])
        #spreadsheet$Design_ID[(i-1)*n.geno + (j-1)*(block.size*n.block.y) + 1:(block.size*n.block.y)] <- as.vector(as.matrix(rep.list[[i]])[,1:n.block.y + (j-1)*(n.block.y)])
        spreadsheet$column[(i-1)*n.geno + (j-1)*(block.size*n.block.y) + 1:(block.size*n.block.y)] <- rep((j-1)*block.size + rep(1:block.size,n.block.y)) + as.numeric(rep.layout=="hor") * (i-1) * n.block.x * block.size
        spreadsheet$row[(i-1)*n.geno + (j-1)*(block.size*n.block.y) + 1:(block.size*n.block.y)]    <- rep(1:n.block.y,each=block.size) + as.numeric(rep.layout=="ver") * (i-1) * n.block.y
      }
    } else {
      spreadsheet$column[(i-1)*n.geno + 1:n.geno] <- rep(1:n.col,n.block.y) + as.numeric(rep.layout=="hor") * (i-1) * n.col
      spreadsheet$row[(i-1)*n.geno + 1:n.geno]    <- rep(1:n.block.y,each=n.col)  + as.numeric(rep.layout=="ver") * (i-1) * n.block.y
    }

    if (rep.layout=="ver") {
      rowNames <- paste('row',1:nrow(layout.list[[i]]) + nrow(layout.list[[i]])*(i-1),sep='')
    } else {
      rowNames <- paste('row',1:nrow(layout.list[[i]]),sep='')
    }



    if (rep.layout=='ver') {
      colNames <- paste('rep',i,'_col_',1:ncol(layout.list[[i]]),sep='')
    } else {
      colNames <- paste('rep',i,'_col_',ncol(layout.list[[1]])*(i-1) + 1:ncol(layout.list[[i]]),sep='')
    }
    row.names(layout.list[[i]]) <- rowNames
    names(layout.list[[i]]) <- colNames
  }

  if (turn.upside.down) {
    for (i in 1:n.rep) {
      layout.list[[i]] <- layout.list[[i]][(nrow(layout.list[[i]])):1,]
    }
    if (rep.layout=="ver") { # change the order of the rep's, if necessary
      layout.list2 <- list()
      for (i in 1:n.rep) {layout.list2[[i]] <- layout.list[[n.rep-i+1]]}
      layout.list <- layout.list2
    }

  }
return(list(rep.list=rep.list,layout.list=layout.list,n.col=n.col,n.geno=n.geno,n.rep=n.rep,rep.layout=rep.layout,spreadsheet=spreadsheet,turn.upside.down=turn.upside.down))
}


MakeDropsDesignFile <- function(wd,ww,fileName,random.order=sample(1:wd$n.geno),geno.list) {
# OBJECTIVE :
# *
# *
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *

  require(XLConnect)

  for (i in 1:wd$n.geno) {
    wd$spreadsheet[wd$spreadsheet$Design_ID==random.order[i],6:8] <- geno.list[i,c(5,3,4)]
  }

  wd$spreadsheet <- data.frame(treatment=rep("WD",nrow(wd$spreadsheet)),wd$spreadsheet)

  wd$layout.list2 <- wd$layout.list

  for (i in 1:wd$n.rep) {
    for (j in 1:ncol(wd$layout.list2[[i]])) {
      wd$layout.list2[[i]][,j] <- as.character(wd$layout.list2[[i]][,j])
    }
  }

  for (i in 1:wd$n.rep) {
    for (j in 1:ncol(wd$layout.list2[[i]])) {
      for (k in 1:nrow(wd$layout.list2[[i]])) {
        wd$layout.list2[[i]][k,j] <- paste(as.character(wd$spreadsheet[min(which(wd$spreadsheet$Design_ID==wd$layout.list[[i]][k,j])),7:9]),collapse="      ")#paste(as.character,as.character,as.character,as.character)
      }
    }
  }

  plot.block.matrix.wd <- matrix("",max(wd$spreadsheet$row),max(wd$spreadsheet$column))
  for (i in 1:nrow(wd$spreadsheet)) {
    plot.block.matrix.wd[wd$spreadsheet$row[i],wd$spreadsheet$column[i]] <- paste(paste('rep',wd$spreadsheet$replicate[i],sep=''),paste('block',wd$spreadsheet$block[i],sep=''),paste('plot',wd$spreadsheet$plot[i],sep=''),sep=',')
  }
  plot.block.matrix.wd <- as.data.frame(plot.block.matrix.wd)
  names(plot.block.matrix.wd) <- paste('column',1:ncol(plot.block.matrix.wd),sep='')
  row.names(plot.block.matrix.wd) <- paste('row',1:nrow(plot.block.matrix.wd),sep='')
  if (wd$turn.upside.down) {plot.block.matrix.wd <- plot.block.matrix.wd[nrow(plot.block.matrix.wd):1,]}

  wd$layout.list3 <- wd$layout.list
  for (i in 1:wd$n.rep) {
    for (j in 1:wd$n.geno) {
      wd$layout.list3[[i]][wd$layout.list[[i]]==random.order[j]] <- geno.list[j,5]
    }
  }
  wd$layout.list <- wd$layout.list3


  #####################################

  # we use the same random.order as for wd !!!!!!!!
  for (i in 1:ww$n.geno) {
    ww$spreadsheet[ww$spreadsheet$Design_ID==random.order[i],6:8] <- geno.list[i,c(5,3,4)]
  }
  ww$spreadsheet <- data.frame(treatment=rep("WW",nrow(ww$spreadsheet)),ww$spreadsheet)

  ww$layout.list2 <- ww$layout.list

  for (i in 1:ww$n.rep) {
    for (j in 1:ncol(ww$layout.list2[[i]])) {
      ww$layout.list2[[i]][,j] <- as.character(ww$layout.list2[[i]][,j])
    }
  }

  for (i in 1:ww$n.rep) {
    for (j in 1:ncol(ww$layout.list2[[i]])) {
      for (k in 1:nrow(ww$layout.list2[[i]])) {
        ww$layout.list2[[i]][k,j] <- paste(as.character(ww$spreadsheet[min(which(ww$spreadsheet$Design_ID==ww$layout.list[[i]][k,j])),7:9]),collapse="        ")#paste(as.character,as.character,as.character,as.character)
      }
    }
  }

  plot.block.matrix.ww <- matrix("",max(ww$spreadsheet$row),max(ww$spreadsheet$column))
  for (i in 1:nrow(ww$spreadsheet)) {
    plot.block.matrix.ww[ww$spreadsheet$row[i],ww$spreadsheet$column[i]] <- paste(paste('rep',ww$spreadsheet$replicate[i],sep=''),paste('block',ww$spreadsheet$block[i],sep=''),paste('plot',ww$spreadsheet$plot[i],sep=''),sep=',')
  }
  plot.block.matrix.ww <- as.data.frame(plot.block.matrix.ww)
  names(plot.block.matrix.ww) <- paste('column',1:ncol(plot.block.matrix.ww),sep='')
  row.names(plot.block.matrix.ww) <- paste('row',1:nrow(plot.block.matrix.ww),sep='')
  if (ww$turn.upside.down) {plot.block.matrix.ww <- plot.block.matrix.ww[nrow(plot.block.matrix.ww):1,]}

  ww$layout.list3 <- ww$layout.list
  for (i in 1:ww$n.rep) {
    for (j in 1:ww$n.geno) {
      ww$layout.list3[[i]][ww$layout.list[[i]]==random.order[j]] <- geno.list[j,5]
    }
  }
  ww$layout.list <- ww$layout.list3

  ###########################################################

  # Load workbook (create if not existing)
  new.sheet <- loadWorkbook(fileName, create = TRUE)
  # Create a worksheet called 'CO2'

  #if (existsSheet(object,name)) {removeSheet(object,sheet)}
  if (existsSheet(object=new.sheet,name="General Info")) {removeSheet(object=new.sheet,sheet="General Info")}
  if (existsSheet(object=new.sheet,name="WW")) {removeSheet(object=new.sheet,sheet="WW")}
  if (existsSheet(object=new.sheet,name="WD")) {removeSheet(object=new.sheet,sheet="WD")}
  if (existsSheet(object=new.sheet,name="WW spreadsheet")) {removeSheet(object=new.sheet,sheet="WW spreadsheet")}
  if (existsSheet(object=new.sheet,name="WD spreadsheet")) {removeSheet(object=new.sheet,sheet="WD spreadsheet")}
  if (existsSheet(object=new.sheet,name="WW full")) {removeSheet(object=new.sheet,sheet="WW full")}
  if (existsSheet(object=new.sheet,name="WD full")) {removeSheet(object=new.sheet,sheet="WD full")}
  if (existsSheet(object=new.sheet,name="WW map")) {removeSheet(object=new.sheet,sheet="WW map")}
  if (existsSheet(object=new.sheet,name="WD map")) {removeSheet(object=new.sheet,sheet="WD map")}

  createSheet(new.sheet, name = "General Info")
  createSheet(new.sheet, name = "WW")
  createSheet(new.sheet, name = "WD")
  createSheet(new.sheet, name = "WW spreadsheet")
  createSheet(new.sheet, name = "WD spreadsheet")
  createSheet(new.sheet, name = "WW full")
  createSheet(new.sheet, name = "WD full")
  createSheet(new.sheet, name = "WW map")
  createSheet(new.sheet, name = "WD map")

  if (ww$rep.layout=="ver") {
    for (i in 1:ww$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WW",data=data.frame(layout=row.names(ww$layout.list[[i]]),ww$layout.list[[i]]),header=T,
                     startCol=1,startRow=1 + (i-1)*(nrow(ww$layout.list[[i]])+2))
    }
  }
  if (ww$rep.layout=="hor") {
    for (i in 1:ww$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WW",data=data.frame(layout=row.names(ww$layout.list[[i]]),ww$layout.list[[i]]),header=T,
                     startRow=1,startCol=1 + (i-1)*(ncol(ww$layout.list[[i]])+2))
    }
  }

  if (wd$rep.layout=="ver") {
    for (i in 1:wd$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WD",data=data.frame(layout=row.names(wd$layout.list[[i]]),wd$layout.list[[i]]),header=T,
                     startCol=1,startRow=1 + (i-1)*(nrow(wd$layout.list[[i]])+2))
    }
  }
  if (wd$rep.layout=="hor") {
    for (i in 1:wd$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WD",data=data.frame(layout=row.names(wd$layout.list[[i]]),wd$layout.list[[i]]),header=T,
                     startRow=1,startCol=1 + (i-1)*(ncol(wd$layout.list[[i]])+2))
    }
  }

  #########

  if (ww$rep.layout=="ver") {
    for (i in 1:ww$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WW full",data=data.frame(layout=row.names(ww$layout.list2[[i]]),ww$layout.list2[[i]]),header=T,
                     startCol=1,startRow=1 + (i-1)*(nrow(ww$layout.list2[[i]])+2))
    }
  }
  if (ww$rep.layout=="hor") {
    for (i in 1:ww$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WW full",data=data.frame(layout=row.names(ww$layout.list2[[i]]),ww$layout.list2[[i]]),header=T,
                     startRow=1,startCol=1 + (i-1)*(ncol(ww$layout.list2[[i]])+2))
    }
  }

  if (wd$rep.layout=="ver") {
    for (i in 1:wd$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WD full",data=data.frame(layout=row.names(wd$layout.list2[[i]]),wd$layout.list2[[i]]),header=T,
                     startCol=1,startRow=1 + (i-1)*(nrow(wd$layout.list2[[i]])+2))
    }
  }
  if (wd$rep.layout=="hor") {
    for (i in 1:wd$n.rep) {
      writeWorksheet(object=new.sheet,sheet="WD full",data=data.frame(layout=row.names(wd$layout.list2[[i]]),wd$layout.list2[[i]]),header=T,
                     startRow=1,startCol=1 + (i-1)*(ncol(wd$layout.list2[[i]])+2))
    }
  }

  writeWorksheet(object=new.sheet,sheet="WD spreadsheet",data=wd$spreadsheet[,-(ncol(wd$spreadsheet))],header=T)
  writeWorksheet(object=new.sheet,sheet="WW spreadsheet",data=ww$spreadsheet[,-(ncol(ww$spreadsheet))],header=T)


  writeWorksheet(object=new.sheet,sheet="WD map",data=data.frame(layout=row.names(plot.block.matrix.wd),plot.block.matrix.wd),header=T)
  writeWorksheet(object=new.sheet,sheet="WW map",data=data.frame(layout=row.names(plot.block.matrix.ww),plot.block.matrix.ww),header=T)


  saveWorkbook(new.sheet)


  return(list(wd=wd,ww=ww,fileName=fileName,random.order=random.order,geno.list=geno.list))
}
#MakeDropsDesignFile(wd=wd,ww=ww,fileName=fileName,random.order=random.order,geno.list=geno.list)

#################################
# partial correlation of y1 and y2, given a matrix X

partial.cor.test <- function(y1,y2,X,qr.decomp=qr(X)) {
  pcor <- cor(reslsfit1(X, y1, qr.decomp),reslsfit1(X, y2, qr.decomp))
  tt   <- sqrt(length(y1) - 3 - ncol(X)) * 0.5 * abs(log((1+pcor)/(1-pcor)))
  pval <- pnorm(tt,lower.tail=T) / 2
  return(c(pcor,tt,pval))
}






#######################################################################################################################################
# STABILITY SELECTION FUNCTIONS


# missing values

# procedures :
# lm
# penalized regression

draw_subsamples <- function(object.names,n.subsample=100,paired=F) {
# OBJECTIVE :
# Drawing pairs of subsamples from ..., of half the sample-size
#
# INPUT :
# * object.names     : the list of genotype-names
# * n.subsample              : total number of subsamples. Should be even; there will be n.subsample/2 pairs of complementary subsamples
# * paired           : do paired subsampling ?
#
# OUTPUT :
# * sub.samples           : a matrix of n.subsample rows, one for each subsample. Every row starts with the observation numbers
#                           i.e. row numbers in pheno.dataframe) of the individuals in the subsample,
#                           then possibly followed by zeros.
# * sub.sample.size       : floor(acc.sample.size/2) (when subsample.acc=T), or
#                           floor(ind.sample.size/2) (when subsample.acc=F)
# * ind.population        : the row numbers (in pheno.dataframe) of the non-missing observations
# * ind.sample.size       : the total number of non-missing observations
# * acc.population        : numbers (i.e. which entry in the vector plant.names) of genotypes
#                           with at least one non-missing observation
# * acc.sample.size       : total number of genotypes with at least one non-missing observation
############################


  if (paired & (round(n.subsample/2)!=n.subsample/2)) {n.subsample  <- n.subsample + 1}

  sample.size <- length(object.names)   # the total number of non-missing observations
  sub.sample.size <- floor(sample.size/2)
  sub.samples         <- matrix(0,n.subsample,sub.sample.size)

  if (paired) {
    for (ss in 1:(n.subsample/2)) {
      sub.samples[2*ss-1,]  <- sort(sample(object.names,size=sub.sample.size))
      sub.samples[2*ss,]    <- sort(sample(setdiff(object.names,sub.samples[2*ss-1,]),size=sub.sample.size))
    }
  } else {
    for (ss in 1:n.subsample) {
      sub.samples[ss,]  <- sort(sample(object.names,size=sub.sample.size))
    }
  }

return(sub.samples)
}

q.elastic.net <- function(y,X,q.select=4,alpha=0.5) {
  require(glmnet)

  regr.obj   <- glmnet(x=as.matrix(X), y=as.numeric(y),family="gaussian",alpha = alpha, standardize = T)
  #regr.obj$beta

  reg.path <- character(0)
  for (lmb in 1:ncol(regr.obj$beta)) {reg.path <- c(reg.path,names(sort(abs(regr.obj$beta[which(regr.obj$beta[,lmb]!=0),lmb]),decreasing=T)))}
  reg.path <- unique(reg.path)
  subsample.selection <- reg.path[1:(min(q.select,length(reg.path)))]

return(subsample.selection)
}


q.elastic.net.full <- function(y,X,alpha=0.5,n.lambda=100) {
  #require(glmnet)
  regr.obj   <- glmnet(x=as.matrix(X), y=as.numeric(y),family="gaussian",alpha = alpha, nlambda = n.lambda, standardize = T)
  regr.obj$beta[regr.obj$beta!=0] <- 1
return(matrix(regr.obj$beta,ncol=ncol(regr.obj$beta)))
}


stability_selection_glmnet <- function(y,X,paired=F,subsamples=draw_subsamples(object.names=rownames(X),paired=paired),Pi=NULL,#procedure='q.elastic.net.full',
                                       ev=1,bt=0,r1 = -1/2,r2 = -1/4,q.select=4,fix.q=F,n.node=1,alpha=0.5,n.lambda=100) {
# add later on : an argument 'procedure' , to generalize
# currently : based on q.elastic.net.full


# OBJECTIVE :

# *
# INPUT :
#
# * X : n x p matrix or data.frame of predictors, without missing values. Must have row- and column names.
# * y : response vector of length n, without missing values
#
# * bt : bound-type (0=Meinshausen & Buhlmann (2010) bound. The other bounds come from Shah & Samworth (2013) : 1=unimodal bound, 2=r-concave bound
#             In case bt=1 or 2, theta = q/p is assumed
#
# * ev : the required bound on the expected number of false positives (V). Ignored if a valid value of Pi is given
#
#
# * Pi : EITHER Pi=NULL, and Pi will be chosen such that E(V) <= ev is satisfied
#        OR  Pi is given a value, and the input ev is ignored (the actual bound for E(V) will be calculated and returned)
#
# * n.node : number of nodes to use. If larger than 1, the parallel package is required (TO DO)
#
# * r1,r2 : values in the r-concavity assumption, needed when bt=2
#
# * q.select : The average number of variables selected for each subsample. If fix.q is TRUE, exactly q.select variables are selected for each subsamples;
#               otherwise the regularization is chosen such that the average number of selected variables is as close as possible to q.select
#
# * fix.q : when TRUE,exactly q.select variables are selected for each subsamples;
#               otherwise the regularization is chosen such that the average number of selected variables is as close as possible to q.select
# *
# * alpha : elastic net parameter
# *
# * n.lambda : size of the grid of the penalty/regularizion parameter
#
# OUTPUT :
# *
# *
# *

  # choice of lambda, q is done in the 'procedure', which is assumed to output a single vector of 0/1 ( a matrix or   list of vectors for several q/lambda values is not (yet) possible)

  #paired=F;subsamples=draw_subsamples(object.names=rownames(X),paired=paired);Pi=NULL;ev=1;bt=0;r1 = -1/2;r2 = -1/4;q.select=4;fix.q=F;n.node=1;alpha=0.5;n.lambda=100

  #subsamples <- draw_subsamples(object.names=rownames(X),n.subsample=100);ev=1;bt=1;paired=T;fix.q=T;Pi=NULL;r1 = -1/2;r2 = -1/4;q.select=4;n.node=1;alpha=0.5;n.lambda=100


  if (class(Pi)!='numeric') {Pi <- NULL}#{class(Pi) <- 'NULL'}
  if ((bt %in% 0:1) & !is.null(Pi)) {stopifnot(Pi > 0.5 & Pi <= 1)} # unless bt==2, Pi needs to be either NULL, or between 0.5 and 1
  if (bt==0) {stopifnot(!paired)}# if bt=0, the subsamples SHOULD NOT BE paired
  if (bt!=0) {stopifnot(paired)} # if bt=1 or 2, the subsamples HAVE TO BE paired

  y <- as.numeric(y)
  X <- as.matrix(X)
  names(y) <- rownames(X)
  p <- ncol(X)
  n.subsample <- nrow(subsamples)

  ###########################

  if (fix.q) {

    selection.matrix <- matrix(0,n.subsample,ncol(X)) # matrix to store the results; subsamples in the rows; vars in the columns; elements will be set to one if a variable was selected for that subsample
    colnames(selection.matrix) <- colnames(X)
    for (i in 1:n.subsample) {
      selection.matrix[i,q.elastic.net(y=y[subsamples[i,]],X=X[subsamples[i,],],q.select=q.select,alpha=alpha)] <- 1
    }
    #q.hat  <- sum(selection.matrix) / n.subsample
    q.opt   <- q.select
    q.hat   <- q.select
    q.index <- 1
    selection.probabilities <- apply(selection.matrix,2,mean)

  } else {

    # array to store the results; subsamples in the 1st dim; variables as 2nd dimension;
    # the grid of regularization parameters as 3rd dimension
    # elements will be set to one if a variable was selected for that subsample
    selection.array  <- array(0,dim=c(n.subsample,ncol(X),n.lambda))
    dimnames(selection.array)[[1]] <- paste('subsample',1:n.subsample,sep='')
    dimnames(selection.array)[[2]] <- colnames(X)
    dimnames(selection.array)[[3]] <- paste('lambda',1:n.lambda,sep='')

    # to do : parallelize
    for (i in 1:n.subsample) {
      selection.array[i,,] <- q.elastic.net.full(y=y[subsamples[i,]],X=X[subsamples[i,],],alpha=alpha,n.lambda=n.lambda)
    }

    q.hat <- apply(selection.array,3,sum) / n.subsample
    q.index <- which.min(abs(q.hat-q.select)) # the index corresponding to the most appropriate regulatization parameter
    q.opt <- as.numeric(q.hat[q.index]) # q.opt : the actual estimate of q
    selection.matrix <- selection.array[,,q.index]
    selection.probabilities <- apply(selection.matrix,2,mean)

  }
  ########################

  if (bt==0) {     # the meinshausen and buhlmann bound
    if (is.null(Pi)) {
      Pi <- 0.5 * (q.opt^2 / (ev*p) + 1)
      if (Pi > 1) {Pi <- 1}
      if (Pi <= 0.5) {Pi <- 0.5001}
    }
    EV <- q.opt^2 / (p * (2*Pi -1))
  }

  if (bt==1) {       # unimodal bound
    if (is.null(Pi)) {
      Pi <- seq(from=0.5 + 2/n.subsample,to=1,by=1/n.subsample)
      EV <- rep(0,length(Pi))
      EV[Pi <= .75] <- as.numeric(.5 * (2*Pi[Pi <= .75]-1-1/n.subsample)^(-1) * q.opt^2 / p)
      EV[Pi >  .75] <- as.numeric(4 * (1-Pi[Pi > .75]+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt^2 / p )
      Pi <- min(Pi[which(EV <= ev)])
      if (length(Pi)==0) {Pi <- 1}
    }
    if (Pi <= .75) {
      EV <- as.numeric(.5 * (2*Pi-1-1/n.subsample)^(-1) * q.opt^2 / p)
    } else {
      EV <- as.numeric(4 * (1-Pi+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt^2 / p )
    }
  }

  if (bt==2) {     # r-concave bounds
    if (is.null(Pi)) {
      Pi <- seq(from=0,to=1,length.out=n.subsample+1)
      EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))
      Pi <- min(Pi[which(EV <= ev)])
      if (length(Pi)==0) {Pi <- 1}
    } else {
      Pi <- round(Pi*n.subsample) / n.subsample # if Pi was specified, it needs to be on the grid 0,1/n.subsample,2/n.subsample,...,1
    }
    EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))[round(Pi*n.subsample)]
  }


return(list(selection=names(which(selection.probabilities>=Pi)),q=q.opt,Pi=Pi,EV=EV,selection.probabilities=sort(selection.probabilities))) # EV[q.index]
}


# test :
# source("D:/willem/Dropbox/research/STATISTICAL_GENETICS/cacao_data_jos/import_cocoa_data.r")
#subsample.matrix <- draw_subsamples(object.names=rownames(X),n.subsample=100)
#(a <- stability_selection_glmnet(y=y,X=X,ev=1,bt=1,paired=T,fix.q=T,subsamples=subsample.matrix))


stability_selection_undirected_graphs <- function(X,paired=F,subsamples=draw_subsamples(object.names=rownames(X),paired=paired),Pi=NULL,#procedure='q.elastic.net.full',
                                                  method='mb',ev=1,bt=0,r1 = -1/2,r2 = -1/4,q.select=20,n.lambda=100) {
# add later on : an argument 'procedure' , to generalize
# currently : based on q.elastic.net.full


# OBJECTIVE :

# *
# INPUT :
#
# * X : n x p matrix or data.frame of predictors, without missing values. Must have row- and column names.
#
# * bt : bound-type (0=Meinshausen & Buhlmann (2010) bound. The other bounds come from Shah & Samworth (2013) : 1=unimodal bound, 2=r-concave bound
#             In case bt=1 or 2, theta = q/p is assumed
#
# * ev : the required bound on the expected number of false positives (V). Ignored if a valid value of Pi is given
#
#
# * Pi : EITHER Pi=NULL, and Pi will be chosen such that E(V) <= ev is satisfied
#        OR  Pi is given a value, and the input ev is ignored (the actual bound for E(V) will be calculated and returned)
#
# * n.node : number of nodes to use. If larger than 1, the parallel package is required (TO DO)
#
# * r1,r2 : values in the r-concavity assumption, needed when bt=2
#
# * q.select : The average number of variables selected for each subsample.
#              The regularization is chosen such that the average number of selected variables is as close as possible to q.select
#
# * method : Graph estimation methods, with 3 options: "mb", "ct" and "glasso". The defaulty value is "mb".
# * n.lambda : the length of the penaly vector, in the huge-package
#
# * alpha : elastic net parameter
# *
# * n.lambda : size of the grid of the penalty/regularizion parameter
#
# OUTPUT :
# *
# *
# *

  # choice of lambda, q is done in the 'procedure', which is assumed to output a single vector of 0/1 ( a matrix or   list of vectors for several q/lambda values is not (yet) possible)

  # paired=T;subsamples=draw_subsamples(object.names=rownames(X),paired=paired);Pi=NULL;method='mb';ev=1;bt=1;r1 = -1/2;r2 = -1/4;q.select=100;n.lambda=100
  # X=f2.matrix;paired=F;subsamples=draw_subsamples(object.names=rownames(X),paired=paired);Pi=NULL;method='mb';ev=1;bt=0;q.select=20;n.lambda=100;r1 = -1/2;r2 = -1/4


  require(huge)

  if (class(Pi)!='numeric') {Pi <- NULL}#{class(Pi) <- 'NULL'}
  if ((bt %in% 0:1) & !is.null(Pi)) {stopifnot(Pi > 0.5 & Pi <= 1)} # unless bt==2, Pi needs to be either NULL, or between 0.5 and 1
  if (bt==0) {stopifnot(!paired)}# if bt=0, the subsamples SHOULD NOT BE paired
  if (bt!=0) {stopifnot(paired)} # if bt=1 or 2, the subsamples HAVE TO BE paired

  X <- as.matrix(X)
  p <- (ncol(X)) * (ncol(X) - 1) / 2 # total number of edges

  n.subsample <- nrow(subsamples)

  huge.obj   <- huge(x=as.matrix(X[subsamples[1,],]), nlambda = n.lambda, method = method)
  path.list  <- huge.obj$path

  for (i in 2:n.subsample) {
    huge.obj   <- huge(x=as.matrix(X[subsamples[i,],]), nlambda = n.lambda, method = method)
    for (j in 1:length(path.list)) {
      path.list[[j]] <- path.list[[j]] + huge.obj$path[[j]]
    }
  }

  for (j in 1:length(path.list)) {
    path.list[[j]] <- path.list[[j]] / n.subsample
  }

  q.hat <- unlist(lapply(path.list,function(M){sum(as.matrix(M)[lower.tri(as.matrix(M))])}))

  q.index <- which.min(abs(q.hat-q.select)) # the index corresponding to the most appropriate regulatization parameter
  q.opt <- as.numeric(q.hat[q.index]) # q.opt : the actual estimate of q

  selection.probabilities <- path.list[[q.index]]

  ########################

  if (bt==0) {     # the meinshausen and buhlmann bound
    if (is.null(Pi)) {
      Pi <- 0.5 * (q.opt^2 / (ev*p) + 1)
      if (Pi > 1) {Pi <- 1}
      if (Pi <= 0.5) {Pi <- 0.5001}
    }
    EV <- q.opt^2 / (p * (2*Pi -1))
  }

  if (bt==1) {       # unimodal bound
    if (is.null(Pi)) {
      Pi <- seq(from=0.5 + 2/n.subsample,to=1,by=1/n.subsample)
      EV <- rep(0,length(Pi))
      EV[Pi <= .75] <- as.numeric(.5 * (2*Pi[Pi <= .75]-1-1/n.subsample)^(-1) * q.opt^2 / p)
      EV[Pi >  .75] <- as.numeric(4 * (1-Pi[Pi > .75]+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt^2 / p )
      Pi <- min(Pi[which(EV <= ev)])
      if (length(Pi)==0) {Pi <- 1}
    }
    if (Pi <= .75) {
      EV <- as.numeric(.5 * (2*Pi-1-1/n.subsample)^(-1) * q.opt^2 / p)
    } else {
      EV <- as.numeric(4 * (1-Pi+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt^2 / p )
    }
  }

  if (bt==2) {     # r-concave bounds
    if (is.null(Pi)) {
      Pi <- seq(from=0,to=1,length.out=n.subsample+1)
      EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))
      Pi <- min(Pi[which(EV <= ev)])
      if (length(Pi)==0) {Pi <- 1}
    } else {
      Pi <- round(Pi*n.subsample) / n.subsample # if Pi was specified, it needs to be on the grid 0,1/n.subsample,2/n.subsample,...,1
    }
    EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))[round(Pi*n.subsample)]
  }

  selection.probabilities <- as.matrix(selection.probabilities)
  colnames(selection.probabilities) <- rownames(selection.probabilities) <- colnames(X)
  selection.probabilities[upper.tri(selection.probabilities)] <- 0
  selection.frame <- which(selection.probabilities>= max(0.2,Pi) ,arr.ind=T)
  selection.frame <- data.frame(variable1=colnames(X)[selection.frame[,1]],variable2=colnames(X)[selection.frame[,2]],
                                selection.probability=selection.probabilities[as.matrix(selection.frame)])
  selection.frame <-  selection.frame[order(selection.frame$selection.probability,decreasing=T),]
  row.names(selection.frame) <- 1:nrow(selection.frame)

  estimated.graph <- selection.probabilities + t(selection.probabilities)
  estimated.graph[estimated.graph >= Pi] <- 1
  estimated.graph <- graph.adjacency(estimated.graph,mode='undirected')

return(list(selection=selection.frame,q=q.opt,Pi=Pi,EV=EV,selection.probabilities=selection.probabilities,estimated.graph=estimated.graph))
}


stability_selection_undirected_graphs_with_two_groups <- function(X,group1,paired=F,subsamples=draw_subsamples(object.names=rownames(X),paired=paired),Pi.inter=NULL,Pi.intra1=NULL,Pi.intra2=NULL,
                                                  method='mb',ev.inter=1,ev.intra1=1,ev.intra2=1,bt=0,r1 = -1/2,r2 = -1/4,q.select.inter=10,q.select.intra1=10,q.select.intra2=10,n.lambda=100) {
# add later on : an argument 'procedure' , to generalize
# currently : based on q.elastic.net.full


# OBJECTIVE :

# *
# INPUT :
#
# * X : n x p matrix or data.frame of predictors, without missing values. Must have row- and column names.
#
# * group1 : column numbers (of X) corresponding to the variables that are in the first group. The others are assumed to be in the second group
#   Accordingly, the edges are split into two groups :
#   - intra-group edges : edges between two variables that are in the same group (either group 1 or group2)
#   - inter-group edges : edges between two variables that are in different groups
#
# * bt : bound-type (0=Meinshausen & Buhlmann (2010) bound. The other bounds come from Shah & Samworth (2013) : 1=unimodal bound, 2=r-concave bound
#             In case bt=1 or 2, theta = q/p is assumed
#
# * ev : the required bound on the expected number of false positives (V). Ignored if a valid value of Pi is given
#
#
# * Pi : EITHER Pi=NULL, and Pi will be chosen such that E(V) <= ev is satisfied
#        OR  Pi is given a value, and the input ev is ignored (the actual bound for E(V) will be calculated and returned)
#
# * n.node : number of nodes to use. If larger than 1, the parallel package is required (TO DO)
#
# * r1,r2 : values in the r-concavity assumption, needed when bt=2
#
# * q.select : The average number of variables selected for each subsample.
#              The regularization is chosen such that the average number of selected variables is as close as possible to q.select
#
# * method : Graph estimation methods, with 3 options: "mb", "ct" and "glasso". The defaulty value is "mb".
# * n.lambda : the length of the penaly vector, in the huge-package
#
# * alpha : elastic net parameter
# *
# * n.lambda : size of the grid of the penalty/regularizion parameter
#
# OUTPUT :
# *
# *
# *

  # choice of lambda, q is done in the 'procedure', which is assumed to output a single vector of 0/1 ( a matrix or   list of vectors for several q/lambda values is not (yet) possible)

  # paired=T;subsamples=draw_subsamples(object.names=rownames(X),paired=paired);Pi=NULL;method='mb';ev=1;bt=1;r1 = -1/2;r2 = -1/4;q.select=100;n.lambda=100
  # group1=sensory.set; X=f2.matrix;paired=F;subsamples=draw_subsamples(object.names=rownames(X),paired=paired);Pi=NULL;method='mb';ev=1;bt=0;q.select=20;n.lambda=100;r1 = -1/2;r2 = -1/4;q.select.inter=10;q.select.intra1=10;q.select.intra2=10;ev.inter=1;ev.intra1=1;ev.intra2=1;Pi.inter=NULL;Pi.intra1=NULL;Pi.intra2=NULL


  require(huge)

  group2 <- setdiff(1:ncol(X),group1)

  if (class(Pi)!='numeric') {Pi <- NULL}#{class(Pi) <- 'NULL'}
  if ((bt %in% 0:1) & !is.null(Pi)) {stopifnot(Pi > 0.5 & Pi <= 1)} # unless bt==2, Pi needs to be either NULL, or between 0.5 and 1

  if (bt==0) {stopifnot(!paired)}# if bt=0, the subsamples SHOULD NOT BE paired
  if (bt!=0) {stopifnot(paired)} # if bt=1 or 2, the subsamples HAVE TO BE paired

  X <- as.matrix(X)
  p <- (ncol(X)) * (ncol(X) - 1) / 2 # total number of edges
  p.inter  <- length(group1) * length(group2)
  p.intra1 <- (length(group1)) * (length(group1) - 1) / 2
  p.intra2 <- (length(group2)) * (length(group2) - 1) / 2


  n.subsample <- nrow(subsamples)

  huge.obj   <- huge(x=as.matrix(X[subsamples[1,],]), nlambda = n.lambda, method = method)
  path.list  <- huge.obj$path

  for (i in 2:n.subsample) {
    huge.obj   <- huge(x=as.matrix(X[subsamples[i,],]), nlambda = n.lambda, method = method)
    for (j in 1:length(path.list)) {
      path.list[[j]] <- path.list[[j]] + huge.obj$path[[j]]
    }
  }

  for (j in 1:length(path.list)) {
    path.list[[j]] <- path.list[[j]] / n.subsample
  }

  path.list.inter <- lapply(path.list,function(M){M[group1,group2]})
  path.list.intra1<- lapply(path.list,function(M){M[group1,group1]})
  path.list.intra2<- lapply(path.list,function(M){M[group2,group2]})

  q.hat <- unlist(lapply(path.list,function(M){sum(as.matrix(M)[lower.tri(as.matrix(M))])}))
  q.hat.inter  <- unlist(lapply(path.list.inter, function(M){sum(as.matrix(M))}))
  q.hat.intra1 <- unlist(lapply(path.list.intra1,function(M){sum(as.matrix(M)[lower.tri(as.matrix(M))])}))
  q.hat.intra2 <- unlist(lapply(path.list.intra2,function(M){sum(as.matrix(M)[lower.tri(as.matrix(M))])}))
  #q.hat - (q.hat.inter + q.hat.intra1 + q.hat.intra2)

  q.index <- which.min(abs(q.hat-q.select)) # the index corresponding to the most appropriate regulatization parameter
  q.index.inter <- which.min(abs(q.hat.inter-q.select.inter)) # the index corresponding to the most appropriate regulatization parameter
  q.index.intra1 <- which.min(abs(q.hat.intra1-q.select.intra1)) # the index corresponding to the most appropriate regulatization parameter
  q.index.intra2 <- which.min(abs(q.hat.intra2-q.select.intra2)) # the index corresponding to the most appropriate regulatization parameter

  q.opt <- as.numeric(q.hat[q.index]) # q.opt : the actual estimate of q
  q.opt.inter  <- as.numeric(q.hat.inter[q.index.inter]) # q.opt : the actual estimate of q
  q.opt.intra1 <- as.numeric(q.hat.intra1[q.index.intra1]) # q.opt : the actual estimate of q
  q.opt.intra2 <- as.numeric(q.hat.intra2[q.index.intra2]) # q.opt : the actual estimate of q


  selection.probabilities <- path.list[[q.index]]
  selection.probabilities.inter <- path.list.inter[[q.index.inter]]
  selection.probabilities.intra1 <- path.list.intra1[[q.index.intra1]]
  selection.probabilities.intra2 <- path.list.intra2[[q.index.intra2]]

  ########################

  if (bt==0) {     # the meinshausen and buhlmann bound

    if (is.null(Pi.inter)) {
      Pi.inter <- 0.5 * (q.opt^2 / (ev.inter*p.inter) + 1)
      if (Pi.inter > 1) {Pi.inter <- 1}
      if (Pi.inter <= 0.5) {Pi.inter <- 0.5001}
    }
    EV.inter <- q.opt.inter^2 / (p.inter * (2*Pi.inter -1))

    if (is.null(Pi.intra1)) {
      Pi.intra1 <- 0.5 * (q.opt.intra1^2 / (ev.intra1*p.intra1) + 1)
      if (Pi.intra1 > 1) {Pi.intra1 <- 1}
      if (Pi.intra1 <= 0.5) {Pi.intra1 <- 0.5001}
    }
    EV.intra1 <- q.opt.intra1^2 / (p.intra1 * (2*Pi.intra1 -1))

    if (is.null(Pi.intra2)) {
      Pi.intra2 <- 0.5 * (q.opt.intra2^2 / (ev.intra2*p.intra2) + 1)
      if (Pi.intra2 > 1) {Pi.intra2 <- 1}
      if (Pi.intra2 <= 0.5) {Pi.intra2 <- 0.5001}
    }
    EV.intra2 <- q.opt.intra2^2 / (p.intra2 * (2*Pi.intra2 -1))

  }

  if (bt==1) {       # unimodal bound

    if (is.null(Pi.inter)) {
      Pi.inter <- seq(from=0.5 + 2/n.subsample,to=1,by=1/n.subsample)
      EV.inter <- rep(0,length(Pi.inter))
      EV.inter[Pi.inter <= .75] <- as.numeric(.5 * (2*Pi.inter[Pi.inter <= .75]-1-1/n.subsample)^(-1) * q.opt.inter^2 / p.inter)
      EV.inter[Pi.inter >  .75] <- as.numeric(4 * (1-Pi.inter[Pi.inter > .75]+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.inter^2 / p.inter )
      Pi.inter <- min(Pi.inter[which(EV.inter <= ev.inter)])
      if (length(Pi.inter)==0) {Pi.inter <- 1}
    }

    if (is.null(Pi.intra1)) {
      Pi.intra1 <- seq(from=0.5 + 2/n.subsample,to=1,by=1/n.subsample)
      EV.intra1 <- rep(0,length(Pi.intra1))
      EV.intra1[Pi.intra1 <= .75] <- as.numeric(.5 * (2*Pi.intra1[Pi.intra1 <= .75]-1-1/n.subsample)^(-1) * q.opt.intra1^2 / p.intra1)
      EV.intra1[Pi.intra1 >  .75] <- as.numeric(4 * (1-Pi.intra1[Pi.intra1 > .75]+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.intra1^2 / p.intra1 )
      Pi.intra1 <- min(Pi.intra1[which(EV.intra1 <= ev.intra1)])
      if (length(Pi.intra1)==0) {Pi.intra1 <- 1}
    }

    if (is.null(Pi.intra2)) {
      Pi.intra2 <- seq(from=0.5 + 2/n.subsample,to=1,by=1/n.subsample)
      EV.intra2 <- rep(0,length(Pi.intra2))
      EV.intra2[Pi.intra2 <= .75] <- as.numeric(.5 * (2*Pi.intra2[Pi.intra2 <= .75]-1-1/n.subsample)^(-1) * q.opt.intra2^2 / p.intra2)
      EV.intra2[Pi.intra2 >  .75] <- as.numeric(4 * (1-Pi.intra2[Pi.intra2 > .75]+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.intra2^2 / p.intra2 )
      Pi.intra2 <- min(Pi.intra2[which(EV.intra2 <= ev.intra2)])
      if (length(Pi.intra2)==0) {Pi.intra2 <- 1}
    }


    if (Pi.inter <= .75) {
      EV.inter <- as.numeric(.5 * (2*Pi.inter-1-1/n.subsample)^(-1) * q.opt.inter^2 / p.inter)
    } else {
      EV.inter <- as.numeric(4 * (1-Pi.inter+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.inter^2 / p.inter )
    }

    if (Pi.intra1 <= .75) {
      EV.intra1 <- as.numeric(.5 * (2*Pi.intra1-1-1/n.subsample)^(-1) * q.opt.intra1^2 / p.intra1)
    } else {
      EV.intra1 <- as.numeric(4 * (1-Pi.intra1+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.intra1^2 / p.intra1 )
    }
    if (Pi.intra2 <= .75) {
      EV.intra2 <- as.numeric(.5 * (2*Pi.intra2-1-1/n.subsample)^(-1) * q.opt.intra2^2 / p.intra2)
    } else {
      EV.intra2 <- as.numeric(4 * (1-Pi.intra2+1/n.subsample)* (1+1/(2*n.subsample))^(-1) * q.opt.intra2^2 / p.intra2 )
    }

  }
  #################################################### to do :
  if (bt==2) {     # r-concave bounds
    if (is.null(Pi)) {
      Pi <- seq(from=0,to=1,length.out=n.subsample+1)
      EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))
      Pi <- min(Pi[which(EV <= ev)])
      if (length(Pi)==0) {Pi <- 1}
    } else {
      Pi <- round(Pi*n.subsample) / n.subsample # if Pi was specified, it needs to be on the grid 0,1/n.subsample,2/n.subsample,...,1
    }
    EV <- p * minD(theta=q.opt/p, n.subsample/2, r = c(r1,r2))[round(Pi*n.subsample)]
  }
  ###########################    ok :


  selection.probabilities.inter <- as.matrix(selection.probabilities.inter)
  colnames(selection.probabilities.inter) <- colnames(X)[group2]
  rownames(selection.probabilities.inter) <- colnames(X)[group1]

  selection.probabilities.intra1 <- as.matrix(selection.probabilities.intra1)
  colnames(selection.probabilities.intra1) <- rownames(selection.probabilities.intra1) <- colnames(X)[group1]
  selection.probabilities.intra1[upper.tri(selection.probabilities.intra1)] <- 0

  selection.probabilities.intra2 <- as.matrix(selection.probabilities.intra2)
  colnames(selection.probabilities.intra2) <- rownames(selection.probabilities.intra2) <- colnames(X)[group2]
  selection.probabilities.intra2[upper.tri(selection.probabilities.intra2)] <- 0

  ####  to do :

  #selection.frame <- which(selection.probabilities>= min(0.2,Pi) ,arr.ind=T)
  #selection.frame <- data.frame(variable1=colnames(X)[selection.frame[,1]],variable2=colnames(X)[selection.frame[,2]],
  #                              selection.probability=selection.probabilities[as.matrix(selection.frame)])
  #selection.frame <-  selection.frame[order(selection.frame$selection.probability,decreasing=T),]
  #row.names(selection.frame) <- 1:nrow(selection.frame)
  #
  #estimated.graph <- selection.probabilities
  #estimated.graph[estimated.graph >= Pi] <- 1
  #estimated.graph <- graph.adjacency(estimated.graph,mode='undirected')

  #####

  selection.frame.intra1 <- which(selection.probabilities.intra1>= min(0.2,Pi.intra1) ,arr.ind=T)
  selection.frame.intra1 <- data.frame(variable1=colnames(X)[group1][selection.frame.intra1[,1]],variable2=colnames(X)[group1][selection.frame.intra1[,2]],
                                selection.probability=selection.probabilities.intra1[as.matrix(selection.frame.intra1)])
  selection.frame.intra1 <-  selection.frame.intra1[order(selection.frame.intra1$selection.probability,decreasing=T),]
  row.names(selection.frame.intra1) <- 1:nrow(selection.frame.intra1)

  selection.frame.intra2 <- which(selection.probabilities.intra2>= min(0.2,Pi.intra2) ,arr.ind=T)
  selection.frame.intra2 <- data.frame(variable1=colnames(X)[group2][selection.frame.intra2[,1]],variable2=colnames(X)[group2][selection.frame.intra2[,2]],
                                selection.probability=selection.probabilities.intra2[as.matrix(selection.frame.intra2)])
  selection.frame.intra2 <-  selection.frame.intra2[order(selection.frame.intra2$selection.probability,decreasing=T),]
  row.names(selection.frame.intra2) <- 1:nrow(selection.frame.intra2)

  selection.frame.inter <- which(selection.probabilities.inter>= min(0.2,Pi.inter) ,arr.ind=T)
  selection.frame.inter <- data.frame(variable1=colnames(X)[group1][selection.frame.inter[,1]],variable2=colnames(X)[group2][selection.frame.inter[,2]],
                                selection.probability=selection.probabilities.inter[as.matrix(selection.frame.inter)])
  selection.frame.inter <-  selection.frame.inter[order(selection.frame.inter$selection.probability,decreasing=T),]
  row.names(selection.frame.inter) <- 1:nrow(selection.frame.inter)


  #estimated.graph.intra1 <- rbind(cbind(selection.probabilities.intra1+t(selection.probabilities.intra1),selection.probabilities.inter),cbind(t(selection.probabilities.inter),selection.probabilities.intra2+t(selection.probabilities.intra2)))#selection.probabilities
  #estimated.graph.intra1[estimated.graph >= Pi] <- 1
  #estimated.graph.intra1 <- graph.adjacency(estimated.graph.intra1,mode='undirected')

  estimated.graph <- rbind(cbind(selection.probabilities.intra1+t(selection.probabilities.intra1),selection.probabilities.inter),cbind(t(selection.probabilities.inter),selection.probabilities.intra2+t(selection.probabilities.intra2)))#selection.probabilities
  estimated.graph[group1,group1][estimated.graph[group1,group1] >= Pi.intra1] <- 1
  estimated.graph[group2,group2][estimated.graph[group2,group2] >= Pi.intra2] <- 1
  estimated.graph[group1,group2][estimated.graph[group1,group2] >= Pi.inter]  <- 1
  estimated.graph[group2,group1][estimated.graph[group2,group1] >= Pi.inter]  <- 1
  estimated.graph[estimated.graph<1] <- 0

  estimated.graph <-   estimated.graph[colnames(X),colnames(X)]
  #which( estimated.graph['mouthfeel_contract',]==1)
  estimated.graph <- graph.adjacency(estimated.graph,mode='undirected')
  #plot(estimated.graph, main="metabolic network", vertex.color="yellow", edge.color="black", vertex.size=8, vertex.frame.color="red", vertex.label=colnames(X))

return(list(highest.selection.probabilities=list(group1=selection.frame.intra1,group2=selection.frame.intra2,between.groups=selection.frame.inter),
            q=list(group1=q.opt.intra1,group2=q.opt.intra2,between.groups=q.opt.inter),
            Pi=list(group1=Pi.intra1,group2=Pi.intra2,between.groups=Pi.inter),
            EV=list(group1=EV.intra1,group2=EV.intra2,between.groups=EV.inter),
            selection.probabilities=list(group1=selection.probabilities.intra1,group2=selection.probabilities.intra2,between.groups=selection.probabilities.inter),
            estimated.graph=estimated.graph))
}




#lod.sc<-lod.scores;lod.sc.rest=lod.scores.rest;gm=gamma;fwe.level=FWE.level;max.d=max.disc;bt=BT;tau=pi
#lod.sc=lod.scores;gm=gamma;max.d=max.disc;bt=bound.type;tau=pi;lod.sc.rest=data.frame(NULL)


StabilitySelectionOptimize <- function(lod.sc,lod.sc.rest=data.frame(NULL),
                                       gm=0,evbound.max=0.05,max.d=100,P=N,bt=2,tau=0,pi.min=0.31, pi.max=0.80) {
# OBJECTIVE :
# * Based on a matrix of lod-scores for subsamples, select markers
#   according to some stability selection bound.
# * Selection is based on the lod-scores themselves, not on ranks (as in StabilitySelectionOptimize2)
# * Bounds of the type E(V) <= EV.level, under P_0 (unimodal bound or r-concave bound)
# WHAT about gamma ?
#
# INPUT :
# * lod.sc is the data-frame with lod-scores for the subsamples for which only a selection of the markers was analyzed
#   N.B. lod.sc should have row-names corresponding to the marker names !!
# * lod.sc.rest is the data-frame with lod-scores for the subsamples for which all markers were analyzed
#   (this data-frame ALSO contains the markers of the selection in lod.sc)
# * bt= bound-type (1=unimodal bound, 2=r-concave bound)
# * max.d = maximal number of discoveries
# * P = number of variables (markers) under consideration
# * tau =  a local version of pi, i.e. the fraction of subsamples for which a variable needs to be selected to be in the final
#        (stability) selection. tau is either a fixed value in (0,1], or tau=0, in which case we will optimize over
#        tau =  pi.min,pi.min+1/ns,pi.min+2/ns,....,pi.max
#      'Optimzing' means : maximizing the number of discoveries under the constraint that E(V) <= evbound.max
# * pi.min, pi.max = the minimal/maximal value for tau
#
# OUTPUT :
# * snp.sel
# * ... impact
# * significance.flag : TRUE if at least one marker is included in the selection; FALSE otherwise
# * thr.b : the optimal LOD-threshold
# * evbound : the bound for E(V) when thr.b is used as LOD-threshold (which may be smaller than evbound.max)
# * selected.markers : the names of the selected markers
# * If significance.flag=FALSE, we set thr.b=Inf, evbound=0 and  selected.markers=NULL
# * pi : if the input tau = 0, the optimal value of tau (WHAT if nothing is found?);
#     if the input tau > 0, this value is the returned pi
#
###################################################################################################
# START OF THE FUNCTION

sign.snps <- TRUE
evbound     <- 0
thr.b       <- Inf

if (bt!= 1 & bt != 2) {bt <- 2}
if (bt==1 & tau==0) {tau=.6}
if (bt==2 & tau!=0) {tau=0} # to do : when bt=2, offer the option to choose a fixed tau

ns          <- ncol(lod.sc)
if (ncol(lod.sc.rest)!=0) {ns.rest <- ncol(lod.sc.rest)} else {ns.rest <- 0}
Nls         <- nrow(lod.sc)

pi.min   <- round(pi.min*ns)/ns  # pi.min and pi.max should be on the grid 1/ns,2/ns,...,1;
pi.max   <- round(pi.max*ns)/ns  # if not round them to the nearest point on this grid
tau.grid <- seq(from=ceiling(ns*pi.min)/ns,to=floor(ns*pi.max)/ns,by=1/ns)

# we now define a matrix EV.table, which will contain information on the max.d best markers, where 'best'
# is determined on the basis the tau*100% quantiles of the lod-scores of the markers. Some (or none) of these max.d markers
# will be eventually selected.
# * the first column will contain the tau*100% quantiles of the lod-scores of the markers, from large to small
# * the second column will contain the value of q for that marker, if the lod-threshold is as in the first column
# * In the third and following columns, we calculate whether, according to some stability selection bound (either bt=1 or 2),
#   a marker is selected for a threshold as in the first column. The corresponding entry in EV[3,..] is 1 if this is the case;
#   otherwise zero. If tau has a fixed value (i.e. when the input tau > 0), EVtable has 3 columns.
#   If tau=0, we optimize over tau in tau.grid. For every value in tau.grid there is a column in EVtable. The first
#   2 columns of EVtable are updated for every different tau in tau.grid
#

if (tau==0) {EV.table    <- matrix(0,max.d,2+length(tau.grid))} else {EV.table    <- matrix(0,max.d,3)}

EV.table[,1]<- sort(apply(lod.scores,1,function(x){sort(x)[round(ns*(1-tau))+1]}),decreasing=TRUE)[1:max.d]   # or use ceiling or floor
# N.B. We require that the proportion of 'good' subsamples is strictly larger than tau, not >= !! (see also the comment below)
# N.B.2. Now we DO assume >= ... after the correction round(ns*(1-tau)) --> round(ns*(1-tau))+1
lod.ecdf    <- ecdf(as.numeric(lod.scores))
EV.table[,2]<- (1-lod.ecdf(EV.table[,1]))*Nls
if (ns.rest>0) {
  lod.ecdf.rest <- ecdf(as.numeric(lod.scores.rest))
  EV.table[,2]<- EV.table[,2] + (1-lod.ecdf.rest(EV.table[,1]))*(P-Nls)
}
# Recall that snp.sel is the vector of of snp numbers whose (whole population) p-value is below alpha=0.05;
# other snps are not considered
comparison.vector   <- 1:max.d
if (gm==0)   {
  comparison.vector   <- rep(1,max.d)
  gm                  <- evbound.max
}


if (bt==1) {
  if (tau <= .75) {
    EV.table[,3]<- as.numeric(.5 * (2*tau-1-1/ns)^(-1) * (EV.table[,2])^2 / (P*comparison.vector) <= gm)
  } else {
    EV.table[,3]<- as.numeric(4 * (1-tau+1/ns)* (1+1/(2*ns))^(-1) * (EV.table[,2])^2 / (P*comparison.vector) <= gm)
  }
}

if (bt==2) {
  lowest.quantile   <-  sort(apply(lod.scores,1,function(x){sort(x)[round(ns*(1-pi.max))+1]}),decreasing=TRUE)[max.d]
  highest.quantile  <-  sort(apply(lod.scores,1,function(x){sort(x)[round(ns*(1-pi.min))+1]}),decreasing=TRUE)[1]

  largest.q         <- (1-lod.ecdf(lowest.quantile))*Nls
  smallest.q        <- (1-lod.ecdf(highest.quantile))*Nls

  if (ns.rest>0) {
    largest.q <- largest.q + (1-lod.ecdf.rest(lowest.quantile))*(P-Nls)
    smallest.q <- smallest.q + (1-lod.ecdf.rest(highest.quantile))*(P-Nls)
  }

  theta.grid        <- seq(from=smallest.q/P,to=largest.q/P,length=100)
  interpolation.matrix<- matrix(0,length(theta.grid),length(tau.grid))
  for (ds in 1:length(theta.grid)) {interpolation.matrix[ds,]<- minD(theta.grid[ds],ns/2)[round(tau.grid*ns)]}

  for (tu in 1:length(tau.grid)) {
    EV.table[,1]<- sort(apply(lod.scores,1,function(x){sort(x)[round(ns*(1-tau.grid[tu]))+1]}),decreasing=TRUE)[1:max.d]   # or use ceiling or floor
    EV.table[,2]<- (1-lod.ecdf(EV.table[,1]))*Nls
    if (ns.rest>0) {EV.table[,2]<- EV.table[,2] + (1-lod.ecdf.rest(EV.table[,1]))*(P-Nls)}
    q.appr  <- approxfun(x=theta.grid,y=interpolation.matrix[,tu])
    EV.table[,2+tu] <- as.numeric(P * q.appr(EV.table[,2]/P) / comparison.vector <= gm) # minD(EV.table[ds,2]/P,ns/2)
  }
}


# If for at least one row (i.e. one value of the threshold) EV.table[,3] is one, then
# set thr.b equal to EV.table[i,1], where i is the largest index such that EV.table[i,3]=1
if (bt==1) {
  if (max(EV.table[,3])>0) {
    if (sum(EV.table[,3])<max.d) {
      best.row<- min(which(EV.table[,3]==0))-1
      thr.b    <- EV.table[best.row,1]
    } else {
      best.row<- max.d
      thr.b  <- EV.table[max.d,1]
    }
  } else {
    sign.snps <- FALSE
  }
}

if (bt==2) {
  if (max(EV.table[,-c(1:2)])>0) {
    d.opt         <- apply(EV.table[,-c(1:2)],2,FUN=function(x){suppressWarnings(min(which(x==0)))-1})
    d.opt[d.opt==Inf] <- max.d
    best.col      <- length(d.opt)-which.max(d.opt[(length(d.opt)):1])-1  # in case of ties, choose the highest possible value of tau...
    tau           <- tau.grid[best.col] #seq(from=1/ns,to=1,length=ns)[best.col]
    # for the chosen tau, recompute exactly...
    EV.table[,1]<- sort(apply(lod.scores,1,function(x){sort(x)[round(ns*(1-tau))+1]}),decreasing=TRUE)[1:max.d]   # or use ceiling or floor
    EV.table[,2]<- (1-lod.ecdf(EV.table[,1]))*Nls
    if (ns.rest>0) {EV.table[,2]<- EV.table[,2] + (1-lod.ecdf.rest(EV.table[,1]))*(P-Nls)}
    minD.vector <- rep(0,max.d)
    for (ds in 1:max.d) {minD.vector[ds]  <- minD(EV.table[ds,2]/P,ns/2)[which(seq(from=1/ns,to=1,length=ns)==tau)]}
    EV.table[,3]  <- as.numeric(P * minD.vector / comparison.vector <= gm) #
    #
    if (sum(EV.table[,3])<max.d) {
      best.row<- min(which(EV.table[,3]==0))-1
      thr.b    <- EV.table[best.row,1]
    } else {
      best.row<- max.d
      thr.b  <- EV.table[max.d,1]
    }
    #if (sum(EV.table[,3])<max.d) {thr.b    <- EV.table[min(which(EV.table[,3]==0))-1,1]} else {thr.b  <- EV.table[max.d,1]}
    #thr.b     <- EV.table[min(which(EV.table[,2+best.col]==0))-1,1]
  } else {
    sign.snps <- FALSE
  }
}

if (sign.snps)   {

  stab.sel        <- apply(lod.sc,1,function(x){as.numeric(sort(x)[round(ns*(1-tau))+1]>=thr.b)})

  if (bt==1) { # can be simplified, using objects computed before ??  e.g. EVtable + best.row
    if (tau <= .75) {
      if (ns.rest>0) {
        (1-lod.ecdf.rest(thr.b))*(P-Nls)
        evbound         <- ((1-lod.ecdf(thr.b))*Nls+(1-lod.ecdf.rest(thr.b))*(P-Nls))^2 / (P*2*(2*tau-1-1/ns))
        rest.impact     <- ((1-lod.ecdf.rest(thr.b))*(P-Nls))/((1-lod.ecdf(thr.b))*Nls+(1-lod.ecdf.rest(thr.b))*(P-Nls))
      } else {
        evbound         <- ((1-lod.ecdf(thr.b))*Nls)^2 / (P*2*(2*tau-1-1/ns))
      }
    } else {
      if (ns.rest>0) {
        evbound         <- 4* (1-tau + 1/ns) * ((1-lod.ecdf(thr.b))*Nls+(1-lod.ecdf.rest(thr.b))*(P-Nls))^2 / (P*(1+1/(2*ns)))
        rest.impact     <- ((1-lod.ecdf.rest(thr.b))*(P-Nls))/((1-lod.ecdf(thr.b))*Nls+(1-lod.ecdf.rest(thr.b))*(P-Nls))
      } else {
        evbound         <- 4* (1-tau + 1/ns) * ((1-lod.ecdf(thr.b))*Nls)^2 / (P*(1+1/(2*ns)))
      }
    }
  }

  if (bt==2) {
    evbound         <-     P * minD(EV.table[best.row,2]/P,ns/2)[best.col]
  }

}

marker.names<- row.names(lod.sc)

if(sign.snps) {
  return(list(significance.flag=sign.snps,evbound=evbound,thr.b,selected.markers=marker.names[stab.sel==1],pi=tau))
} else {
  return(list(significance.flag=sign.snps,evbound=evbound,thr.b,selected.markers=NULL,pi=tau))
}

}
######### end of function StabilitySelectionOptimize

StabilitySelectionOptimize2 <- function(lod.sc=data.frame(NULL),lod.sc.rest=data.frame(NULL),evbound.max=0.05,
                                        max.d=100,P=N,tau=0,lod.sel=NULL,pi.min=0.50) {
# Objective : ??
#
# lod.sc is the data-frame with lod-scores for the subsamples for which only a selection of the markers was analyzed
# lod.sc.rest is the data-frame with lod-scores for the subsamples for which all markers were analyzed
#   (this data-frame ALSO contains the markers of the selection in lod.sc)
# max.d = maximal number of discoveries
# P = number of variables (markers) under consideration
# tau =  pi
# lod.sel =
# pi.min =

sign.snps <- TRUE
evbound     <- 0
thr.b       <- 0
if (ncol(lod.sc.rest)!=0) {
  lod.sc <- lod.sc.rest
  rm(lod.sc.rest)
  gc()
}
ns          <- ncol(lod.sc)
D.name      <- paste("D.table.",P,".",max.d,".",ns,".RData",sep="")

if (file.exists(D.name)) {
    load(file=D.name)
} else {
    D.matrix  <- matrix(0,max.d,ns)
    for (ds in 1:max.d) {D.matrix[ds,]  <- minD(ds/P,ns/2)}
    save(D.matrix,file=D.name)
}

########
#tau.min <- 0.31; tau.max <- 0.80

if (tau==0) {
  #pi.min  <- 1/ns; pi.max <- 1
  tau.grid <- seq(from=ceiling(ns*pi.min)/ns,to=floor(ns*pi.max)/ns,by=1/ns)
  EV.table <- ((D.matrix * P) <= evbound.max)  # matrix(0,max.d,length(tau.grid))
  d.opt         <- apply(EV.table,2,FUN=function(x){suppressWarnings(min(which(x==FALSE)))-1})
  d.opt[d.opt==Inf] <- max.d
  ev.opt  <- rep(0, length(d.opt))
  #
  # d.opt is a vector containing (for all values in tau.grid)
  # the number of snps that can be taken into account under the constraint that E(V) <= evbound.max (i.e. we look at the ... best snps)
  # ev.opt is a vector containing (for all values in tau.grid)
  # the number of snps that are actually selected when looking at the d.opt[...] best snps in every sub.sample, for tau=tau.grid[...]
  #lod.sc2 <- lod.sc
  #lod.sc2[lod.sc2>max.disc]<- NA
  #lod.sc2 <- lod.sc2[apply(lod.sc2,1,FUN=function(x){sum(!is.na(x))})!=0,]
  #lod.sc2[is.na(lod.sc2)]<-max.disc+1
  ##lod.sc[lod.sc>max.disc]<- NA
  ###lod.selection <- (apply(lod.sc,1,FUN=function(x){sum(!is.na(x))})!=0)
  ###lod.sc <- lod.sc[lod.selection,]
  #
  lod.sc[is.na(lod.sc)]<-max.disc+1
  #
  #for (tu in 1:length(tau.grid)) {if (d.opt[tu]>0) {ev.opt[tu]  <- sum(apply(lod.sc2,1,FUN=function(x){mean(x<=d.opt[tu],na.rm=T)}) >= tau.grid[tu])}}
  for (tu in 1:length(tau.grid)) {
    if (d.opt[tu]>0) {ev.opt[tu]  <- sum(apply(lod.sc,1,FUN=function(x){mean(x<=d.opt[tu],na.rm=T)}) >= tau.grid[tu])}
  }
  # a temporary solution : require that tau is at least pi.min
  ev.opt[1:(round(ns*pi.min))] <- 0
  #if(ns<150) {ev.opt[ev.opt>0][d.opt[ev.opt>0]/ev.opt[ev.opt>0]>100] <- 0} else  {ev.opt[ev.opt>0][d.opt[ev.opt>0]/ev.opt[ev.opt>0]>20] <- 0}

  if (sum(ev.opt>0)) {
    best.col      <- length(ev.opt)-which.max(ev.opt[(length(ev.opt)):1])+1  # in case of ties, choose the highest possible value of tau...
    tau           <- tau.grid[best.col] #seq(from=1/ns,to=1,length=ns)[best.col]
    evbound       <- ev.opt[best.col]
    thr.b         <- d.opt[best.col]
    #stab.sel      <- as.numeric(apply(lod.sc,1,FUN=function(x){mean(x<=thr.b) >= tau})) # NB not length P..., unless NS-NS.rest==0
    stab.sel      <- rep(0,P)
    stab.sel[lod.sel] <- as.numeric(apply(lod.sc,1,FUN=function(x){mean(x<=thr.b) >= tau})) # NB not length P..., unless NS-NS.rest==0
    sign.snps <- TRUE
  } else {
    sign.snps <- FALSE
  }
} else {
  EV.table  <- ((D.matrix * P) <= evbound.max)
  tau       <- round(ns*tau)/ns
  best.col  <- round(ns*tau)
  d.opt     <- apply(EV.table,2,FUN=function(x){suppressWarnings(min(which(x==FALSE)))-1})[best.col]
  if (d.opt==Inf) {d.opt <- max.d}
  ev.opt    <- 0
  if (d.opt>0) {ev.opt  <- sum(apply(lod.sc,1,FUN=function(x){mean(x<=d.opt)}) >= tau)}

  if (ev.opt>0) {
    evbound       <- ev.opt
    thr.b         <- d.opt
    stab.sel      <- as.numeric(apply(lod.sc,1,FUN=function(x){mean(x<=thr.b) >= tau})) # NB not length P...
    sign.snps <- TRUE
  } else {
      sign.snps <- FALSE
  }
}

marker.names  <- row.names(lod.sc)
if(sign.snps) {
  return(list(significance.flag=sign.snps,evbound=evbound,thr.b,selected.markers=marker.names[stab.sel==1],pi=tau))
} else {
  return(list(significance.flag=sign.snps,evbound=evbound,thr.b,selected.markers=NULL,pi=tau))
}

}
######### end of function StabilitySelectionOptimize

######
# CODE BY RAJEN SHAH
r.TailProbs <- function(eta, B, r) {
# TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
# We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
# s is -1/r
  	MAXa <- 100000
	MINa <- 0.00001
	s <- -1/r
	etaB <- eta * B
	k.start <- (ceiling(2 * etaB) + 1)
	if(k.start > B) stop("eta is too large")

	a.vec <- rep(MAXa,B)

	Find.a <- function(prev.a) uniroot(Calc.a, lower = MINa, upper = prev.a, tol = .Machine$double.eps^0.75)$root
	Calc.a <- function(a) {
		denom <- sum((a + 0:k)^(-s))
		num <- sum((0:k) * (a + 0:k)^(-s))
		num / denom - etaB
	}

	for(k in k.start:B) a.vec[k] <- Find.a(a.vec[k-1])

	OptimInt <- function(a) {
		num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
		denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
		1 - num / denom
	}

	output <- rep(1, B)

	prev.k <- k.start
	for(t in k.start:B) {
		cur.optim <- rep(0, B)
		for (k in prev.k:(B-1)) cur.optim[k] <- optimize(f=OptimInt, lower = a.vec[k+1], upper = a.vec[k], maximum  = TRUE)$objective
		output[t] <- max(cur.optim)
		prev.k <- which.max(cur.optim)
	}
	return(output)
}

minD <- function(theta, B, r = c(-1/2, -1/4)) {
  pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])), r.TailProbs(theta, 2*B, r[2]))
}

# END of code by Rajen Shah


# END OF STABILITY SELECTION   SECTION
###########################################################################################################################################
# Now : the obsolete emma package

emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == length(snps))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}

emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}

emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}

emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)

  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                complete=TRUE)[,c(1:(t-q),(t+1):n)]))
}

emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(lambda+delta))))-sum(log(xi+delta))) )
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  delta <- exp(logdelta)
#  stopifnot(length(lambda) == length(etas.1))
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(xi.1+delta))+(n-t)*logdelta)) )
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(n*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/(xi+delta))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(n*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/(xi.1+delta))+(n-t)/delta) ) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(lambda+delta))+(n-t)*logdelta)) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(nq*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/ldelta)+(n-t)/delta)) )
}

emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)


    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))
    dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim

    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,t,m) + matrix(delta,t,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    #LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)+etas.2.sq/delta))-colSums(log(Xis))+(n-t)*log(deltas))
    dLL <- 0.5*delta*(n*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Xis)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/n
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/n
  }
  maxve <- maxva*maxdelta

  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.MLE.noX <- function(y, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL)
{
  n <- length(y)
  t <- nrow(K)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)

  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    etas <- crossprod(eig.L$vectors,y)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n,m)
    LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Xis)))-colSums(log(Xis)))
    dLL <- 0.5*delta*(n*colSums(Etasq/(Xis*Xis))/colSums(Etasq/Xis)-colSums(1/Xis))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    #print(dLL)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.L$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.L$values,etas,eig.L$values))
    }

    for( i in 1:(m-1) )
      {
        #if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.L$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.L$values, etas, eig.L$values))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    etas <- crossprod(eig.L$vectors,y)
    etas.1 <- etas[1:t]
    etas.2 <- etas[(t+1):n]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim

    m <- length(logdelta)
    delta <- exp(logdelta)
    Xis <- matrix(eig.L$values,t,m) + matrix(delta,t,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t,m)
    #LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)+etas.2.sq/delta))-colSums(log(Xis))+(n-t)*log(deltas))
    dLL <- 0.5*delta*(n*(colSums(Etasq/(Xis*Xis))+etas.2.sq/(delta*delta))/(colSums(Etasq/Xis)+etas.2.sq/delta)-(colSums(1/Xis)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.L$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.L$values,etas.1,eig.L$values,n,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.L$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.L$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.L$values+maxdelta))/n
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.L$values+maxdelta))+etas.2.sq/maxdelta)/n
  }
  maxve <- maxva*maxdelta

  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL) {
#y=as.numeric(pheno.frame[,5]);X=cbind(as.matrix(rep(1,n.available*n.rep)),pheno.frame[,-(1:5)]); K=K[pheno.frame$genotype,pheno.frame$genotype]
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  X <- as.matrix(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) { #   min(svd(X)$d)
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta

  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.ML.LRT <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }
  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }

  g <- nrow(ys)
  n <- ncol(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  q1 <- q0 + 1

  if ( !ponly ) {
    ML1s <- matrix(nrow=m,ncol=g)
    ML0s <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)
  ML0 <- vector(length=g)

  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    for(i in 1:g) {
      ML0[i] <- emma.MLE(ys[i,],X0,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R0)$ML
    }

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
          ML0s[i,] <- ML0s[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0[vids,,drop=FALSE],xs[i,vids])
          eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          nr <- sum(vrows)
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)
        }

        for(j in 1:g) {
          if ( nv == t ) {
            MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
#            MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
            if (!ponly) {
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])

          }
          else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              MLE0 <- emma.MLE(ys[j,vids],X0[vids,,drop=FALSE],K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vids],X,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
            }
            else {
              if ( nr == n ) {
                MLE1 <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L)
              }
              else {
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids],K[vids,vids])
                MLE0 <- emma.MLE(ys[j,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
                MLE1 <- emma.MLE(ys[j,vrows],X,K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
              }
            }
            if (!ponly) {
              ML1s[i,j] <- MLE1$ML
              ML0s[i,j] <- MLE0$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }
        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    for(i in 1:g) {
      vrows <- !is.na(ys[i,])
      if ( is.null(Z) ) {
        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)$ML
      }
      else {
        vids <- colSums(Z[vrows,]>0)

        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp)$ML
      }
    }

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[,i] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
#          print(j)
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,NULL,ngrids,llim,ulim,esp,eig.L,eig.R1)
              }
              else {
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)
              }
            }
            else {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vtids,vtids],Z[vrows,vtids],ngrids,llim,ulim,esp)
              }
            }

            if (!ponly) {
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L(NULL,K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtids],X0[vtids,,drop=FALSE],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtids],X[vtids,],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids]))
              eig.L0 <- emma.eigen.L(Z[vtrows,vtids],K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtrows],X0[vtrows,,drop=FALSE],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtrows],X[vtrows,],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
            }
            if (!ponly) {
              ML1s[i,j] <- MLE1$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
              ML0s[i,j] <- MLE0$ML
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }
        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
  }
}

emma.test <- function(ys, xs, K, Z=NULL, x0s = NULL, X0 = NULL, dfxs = 1, dfx0s = 1, use.MLE = FALSE, use.LRT = FALSE, ngrids = 100, llim = -10, ulim = 10, esp=1e-10, ponly = FALSE)
{
  stopifnot (dfxs > 0)

  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }

  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  nx <- nrow(xs)/dfxs

  if ( is.null(x0s) ) {
    dfx0s = 0
    x0s <- matrix(NA,0,ncol(xs))
  }
  # X0 automatically contains intercept. If no intercept is to be used,
  #    X0 should be matrix(nrow=ncol(ys),ncol=0)
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }

  stopifnot(Z == NULL) # The case where Z is not null is not implemented

  ny <- nrow(ys)
  iy <- ncol(ys)
  ix <- ncol(xs)

  stopifnot(nrow(K) == ix)
  stopifnot(ncol(K) == ix)
  stopifnot(nrow(X0) == iy)

  if ( !ponly ) {
    LLs <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  dfs <- matrix(nrow=m,ncol=g)
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)

  # The case with no missing phenotypes
  if ( sum(is.na(ys)) == 0 ) {
    if ( ( use.MLE ) || ( !use.LRT ) ) {
      eig.L0 <- emma.eigen.L(Z,K)
    }
    if ( dfx0s == 0 ) {
      eig.R0 <- emma.eigen.R(Z,K,X0)
    }
    x.prev <- NULL

    for(i in 1:ix) {
      x1 <- t(xs[(dfxs*(i-1)+1):(dfxs*i),,drop=FALSE])
      if ( dfxs0 == 0 ) {
        x0 <- X0
      }
      else {
        x0 <- cbind(t(x0s[(dfx0s*(i-1)+1):(dfx0s*i),,drop=FALSE]),X0)
      }
      x <-  cbind(x1,x0)
      xvids <- rowSums(is.na(x) == 0)
      nxv <- sum(xvids)
      xv <- x[xvids,,drop=FALSE]
      Kv <- K[xvids,xvids,drop=FALSE]
      yv <- ys[j,xvids]

      if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          dfs[i,] <- dfs[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          stats[i,] <- stats[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        eig.R1 = emma.eigen.R.wo.Z(Kv,xv)

        for(j in 1:iy) {
          if ( ( use.MLE ) || ( !use.LRT ) ) {
            if ( nxv < t ) {
              # NOTE: this complexity can be improved by avoiding eigen computation for identical missing patterns
              eig.L0v <- emma.eigen.L.wo.Z(Kv)
            }
            else {
              eig.L0v <- eig.L0
            }
          }

          if ( use.MLE ) {
            MLE <- emma.REMLE(yv,xv,Kv,NULL,ngrids,llim,ulim,esp,eig.R1)
            stop("Not implemented yet")
          }
          else {
            REMLE <- emma.REMLE(yv,xv,Kv,NULL,ngrids,llim,ulim,esp,eig.R1)
            if ( use.LRT ) {
              stop("Not implemented yet")
            }
            else {
              U <- eig.L0v$vectors * matrix(sqrt(1/(eig.L0v$values+REMLE$delta)),t,t,byrow=TRUE)
              dfs[i,j] <- length(eig.R1$values)
              yt <- crossprod(U,yv)
              xt <- crossprod(U,xv)
              ixx <- solve(crossprod(xt,xt))
              beta <- ixx%*%crossprod(xt,yt)
              if ( dfxs == 1 ) {
                stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
              }
              else {
                model.m <- c(rep(1,dfxs),rep(0,ncol(xv)-dfxs))
                stats[i,j] <-
                  crossprod(crossprod(solve(crossprod(crossprod(iXX,model.m),
                                                      model.m)),
                                      model.m*beta),model.m*beta)

              }
              if ( !ponly ) {
                vgs[i,j] <- REMLE$vg
                ves[i,j] <- REMLE$ve
                REMLs[i,j] <- REMLE$REML
              }
            }
          }
        }
        if ( dfxs == 1 ) {
          ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
        }
        else {
          ps[i,] <- pf(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
        }
      }
    }
  }
  # The case with missing genotypes - not implemented yet
  else {
    stop("Not implemented yet")
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          dfs[i,] <- dfs[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids,drop=FALSE]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            yv <- ys[j,vrows]
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,NULL,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),n,n,byrow=TRUE)
              }
              else {
                eig.L0 <- emma.eigen.L.wo.Z(K[vrows,vrows,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vrows,vrows,drop=FALSE],NULL,ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }
            else {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,Z,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }

            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X[vrows,,drop=FALSE])
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L.wo.Z(K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtids]
              nr <- sum(vtids)
              REMLE <- emma.REMLE(yv,X[vtids,,drop=FALSE],K[vtids,vtids,drop=FALSE],NULL,ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtids,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,,drop=FALSE])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids,drop=FALSE]))
              eig.L0 <- emma.eigen.L.w.Z(Z[vtrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtrows]
              nr <- sum(vtrows)
              REMLE <- emma.REMLE(yv,X[vtrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vtrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtrows,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            yt <- crossprod(U,yv)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)

          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,REMLs=REMLs,stats=stats,dfs=dfs,vgs=vgs,ves=ves))
  }
}

emma.REML.t <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }
  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }

  g <- nrow(ys)
  n <- ncol(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  q1 <- q0 + 1

  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( !ponly ) {
    REMLs <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  dfs <- matrix(nrow=m,ncol=g)
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if ( !ponly ) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          stats[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)

      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          dfs[i,] <- dfs[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          stats[i,] <- stats[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0[vids,,drop=FALSE],xs[i,vids])
          eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)
        }

        for(j in 1:g) {
          if ( nv == t ) {
            REMLE <- emma.REMLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.R1)
            if ( is.null(Z) ) {
              U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),t,t,byrow=TRUE)
              dfs[i,j] <- nv - q1
            }
            else {
              U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              dfs[i,j] <- n - q1
            }
            yt <- crossprod(U,ys[j,])
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)

            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              nr <- sum(vids)
              yv <- ys[j,vids]
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],NULL,ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            else {
              eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids,drop=FALSE],K[vids,vids])
              yv <- ys[j,vrows]
              nr <- sum(vrows)
              tv <- sum(vids)
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],Z[vrows,vids,drop=FALSE],ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-tv)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if (!ponly) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          dfs[i,] <- dfs[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids,drop=FALSE]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            yv <- ys[j,vrows]
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,NULL,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),n,n,byrow=TRUE)
              }
              else {
                eig.L0 <- emma.eigen.L.wo.Z(K[vrows,vrows,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vrows,vrows,drop=FALSE],NULL,ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }
            else {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,Z,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }

            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X[vrows,,drop=FALSE])
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L.wo.Z(K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtids]
              nr <- sum(vtids)
              REMLE <- emma.REMLE(yv,X[vtids,,drop=FALSE],K[vtids,vtids,drop=FALSE],NULL,ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtids,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,,drop=FALSE])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids,drop=FALSE]))
              eig.L0 <- emma.eigen.L.w.Z(Z[vtrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtrows]
              nr <- sum(vtrows)
              REMLE <- emma.REMLE(yv,X[vtrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vtrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtrows,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            yt <- crossprod(U,yv)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)

          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,REMLs=REMLs,stats=stats,dfs=dfs,vgs=vgs,ves=ves))
  }
}


# OBJECTIVE :
# *
# *
# INPUT :
# *
# *
# *
# OUTPUT :
# *
# *
# *

