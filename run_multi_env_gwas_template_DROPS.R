################ TO DO :
#
# MAX.DIAG SHOULD DEPEND ON THE SCALE OF THE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!
#
#

# remove all objects
rm(list=ls())
gc()

#######################################################################
# SETTINGS      (some other settings are defined when calling the functions
#                define_pheno_data and EM_function_FA below)
#######################################################################

##########################################################
# Tasks to perform


run.em   <- T  # if F, read the variance components from the file varcomp.file in results.folder

run.gwas <- T  # if F, read the gwas-results from the file gwas.file in results.folder

# Extra string attached to all output file-names; can also be ''
suffix <- ''

# the following option is still under construction; leave to zero
LOD.thr <- 0   # if larger than zero, it is assumed a GWAS was done previously with the same name
               # .. and GWAS is now only (re)run for markers with -log(p) larger than LOD.thr

##########################################################
# Kinship

reduce.K <- F

n.pca <- 25

##########################################################
# Directories

# here are the phenotypic csv files
data.folder    <- './example_data_DROPS/'

# here are the R-scripts
script.folder  <- "./"

# Results will be written to this folder (needs to exist!)
# should have an / at the end
results.folder <- './example_data_DROPS/'

##########################################################
# Phenotypic data

# Vector of phenotypic files to be analyzed (e.g. one per trait)
#pheno.files <- c('grainnumber_multi_env_gwas_dec2015.csv','anthesis_multi_env_gwas_dec2015.csv')
#
# Each file needs to be a csv file (comma-separated; not with semi-colons !),
# and needs to contain genotypic means for a number of environments or traits
# The first column(s) (there may be more than 1) should contain genotype labels.
# In case of multi-trait data, each file should contain the set of traits that should
# be analyzed together.
# In case of multi-environment data, each file should contain all environments for one trait.
# The scripts are not really made for multi-trait/environment simultaneously, although...
pheno.files <- paste0(data.folder,'anthesis_multi_env_gwas_dec2015.csv')

# In the phenotypic files, do all columns start with the trait names ? (as in DROPS)
read.trait.names.from.file <- T

# If read.trait.names.from.file==F :
# Should correspond to pheno.files
# Is only used to create the names of the output files
# Can be a vector of '', e.g. c('','','') in case pheno.files has length 3
# Only important in case you do a multi-environment GWAS for several traits separately,
#  pheno.files containing the files for each trait
trait.names <- ''

# The first n.omit columns: these are all columns not containing phenotypic data,
# in particular genotype labels
n.omit <- 1

# The column-name containing the genotype/accession/hybrid etc labels.
#   Should be one of the first n.omit (see above) columns
#   Should match with the names in GWAS.obj$plant.names , contained in r.image (see above)
#   Always specify; cannot be omitted if n.omit == 1
geno.ID <- 'Accession_ID'

# select a subset of situations; to select all situations: 0
# If not 0, trait.selection should be a vector of integers, referring to situations
# in the order in which they occur in the file(s), e.g. to select the first 10, set trait.selection <- 1:10
trait.selection <- 0

# Standardize all phenotypic columns ?
standardize <- F

# Should imputation be performed ? Can only be F in case of balanced trials
no.imputation <- TRUE

##########################################################
# Genotypic data

# genotypic data (give full path; drops.RData = 50k ; otherwise drops_600k.RData)
r.image <- './example_data_DROPS/drops_600k.RData'

MAF <- 0.05

take.a.subset.of.markers <- T

# If take.a.subset.of.markers==TRUE, the following is used:
marker.subset <- 94001:96000#sort(sample(1:330000,5000)) #sort(c(94831,sample(1:214051,100)))

# name(s) of SNP(s) that are to be included as covariate (vector of strings)
# default :
snp.covariates.list <- list(); for (i in 1:length(pheno.files)) {snp.covariates.list[[i]] <- ''}

#snp.covariates.list <- list('AX-90548584')

# COVARIATES:
#   THE SAME COVARIATES ARE USED FOR ALL files in pheno.files  !
# file with other covariates
covariate.file <- ''   # covariates_test_file.csv


##########################################################
# EM-parameters (fitting the variance components)

# 1 = unstructured for both Vg and Ve (as in Zhou and Stephens)
# 2 = unstructured for both Vg and Ve (pairwise, as in Furlotte and Eskin)
# 3 = factor-analytic for both Vg and Ve
# 4 = read from file. No snp.covariates allowed !
# 5 (NEW here) : approximate ...
cov.model <- 3

# If cov.model==5:
use.heritability.package <- TRUE

# In case cov.model==2 : should there be environmental correlations ?
# If traits are measured on the same plants, put to FALSE
vE.diag <- T#FALSE

# In case cov.model==4, the name of the file
# Should be and .RData file with matrices called Vg and Ve,
# which should have rownames and colnames corresponding to the trait/environment names
# used in the phenotypic files. It can contain a larger number of traits/environments
# than is analyzed here.
# Should be in data.folder
cov.file.name <- '.....RData'

# In case cov.model==3
# tolerance in the EM-algorithm; it stops when the difference in conditional LL
#   in 2 consecutive iterations drops below tol.em
tol.em <- 0.000001

# In case cov.model==3
# Maximum number of iterations in the EM-algorithm
max.iter.em <- 200000

# In case cov.model==3
# Maximal value of the diagonal elements in the precision matrices Cm and Dm
# (ignoring the low-rank part W W^t). recommended value: 100
max.diag <- 10000

# In case cov.model==3
# Factor analytic models:
# Order of genetic part
m.G <- 1

# In case cov.model==3
# Order of environmental part
m.E <- 1

# In case cov.model==3
# Should there be an extra diagonal part in the model for Cm ? (DEFAULT: TRUE)
Cm.het <- TRUE

# In case cov.model==3
# Should there be an extra diagonal part in the model for Dm ? (DEFAULT: TRUE)
Dm.het <- TRUE

# In case cov.model==3
stopifdecreasing <- TRUE

# In case cov.model==3
compute.log.lik <- TRUE

########################################################################################
# START OF THE PROGRAM
########################################################################################

library(MASS)     # for the ginv function
library(mice)
library(corpcor)
library(Matrix)
library(qqman)

##############################################################################################
#setwd(script.folder)

#source(file="EM_functions.R")
#source(file="general_functions.R")
#source(file="newton_raphson_functions_diag.R")

#source(file="newton_raphson_functions.R")

#source(file="updateC.R")
#source(file="EM_function_with_indices.R")

#source(file="update_LS.R")
#source(file="EM_function_LS.R")

#source(file="update_FA.R")
#source(file="update_FA_homogeneous_var.R")
#source(file="EM_function_FA.R")

#source(file='gwas_functions.R')

#source(file="read_and_impute_pheno_file_FUNCTION.R")

#source(file="TransformGwasObject.R")

#source(file="asreml_unstructured_pairwise.R")

#setwd('general_functions/')

#for (file.name in list.files()) {
#  source(file = file.name)
#}
######################################################################

if (cov.model == 5) {
  #stopifnot(length(snp.covariates.list)==0)
  stopifnot(covariate.file=='')
  if (use.heritability.package) {library(heritability)}
}

#########################################################################
#pheno.file = pheno.files[1]

for (pheno.file in pheno.files) {
#pheno.file = pheno.files[1]

  snp.covariates <- snp.covariates.list[[which(pheno.file==pheno.files)]]

  if (cov.model==2) {stopifnot(snp.covariates=='')}

  if (cov.model==4) {stopifnot(snp.covariates=='')}

  load(r.image)

  ###

  if (take.a.subset.of.markers==TRUE) {

    GWAS.obj.original <- GWAS.obj

    #marker.subset

    if (snp.covariates!='') {
      if (length(which(rownames(GWAS.obj$markers) %in% snp.covariates)) > 0) {
        marker.subset <- sort(c(marker.subset, which(rownames(GWAS.obj$markers) %in% snp.covariates)))
        cat('Co-factors have been added to the marker-subset \n')
      }
    }

    GWAS.obj <- TransformGwasObject(gwas.obj=GWAS.obj,
                                    marker.selection=marker.subset,
                                    modify.external=F)

  }

  ###

  #setwd(data.folder)

  if (cov.model==5) {
    standardize <- FALSE
  }

  pheno <- define_pheno_data(pheno.file=pheno.file,
                                standardize=standardize,
                                covariate.file=covariate.file,
                                snp.covariates=snp.covariates,
                                r.image=r.image,
                                n.impute=11,no.imputation=no.imputation,
                                trait.selection=trait.selection,
                                n.omit=n.omit,
                                geno.ID=geno.ID,
                                read.trait.name.from.file=read.trait.names.from.file,
                                trait.name=trait.names[which(pheno.files==pheno.file)]
                                )

  Y <- pheno$Y

  if (cov.model==5) {
    for (j in 1:ncol(Y)) {
      Y[,j] <- Y[,j] - mean(Y[,j])
    }
  }

  Y.unstandardized <- pheno$Y.unstandardized

  K <- pheno$K

  X <- pheno$X

  # run lines in new_stuff.R  !

  if (snp.covariates[1]!='') {
    if ((ncol(X)-length(snp.covariates))==0) {
      X.red <- matrix(0,nrow(X),0)
      rownames(X.red) <- rownames(X)
    } else {
      X.red <- as.matrix(X[,1:(ncol(X)-length(snp.covariates))])
    }
  }

  trait.name <- pheno$trait.name

  #setwd(results.folder)

  write.csv(Y.unstandardized,file=paste0(results.folder,trait.name,suffix,'_imputed_unstandardized.csv'),quote=F)

  varcomp.file <- paste0('varcomp_',trait.name,'_',m.G,'_',m.E,suffix,'.RData')


  if (cov.model==3) {
    gwas.file <- paste0(results.folder,'gwas_',trait.name,'_',m.G,'_',m.E,suffix,'.RData')
    gwas.csv  <- paste0(results.folder,'gwas_effects_',trait.name,'_',m.G,'_',m.E,suffix,'.csv')
    gwas.csv.t<- paste0(results.folder,'gwas_Tstat_',trait.name,'_',m.G,'_',m.E,suffix,'.csv')
  } else {
    gwas.file <- paste0(results.folder,'gwas_',trait.name,'_',suffix,'.RData')
    gwas.csv  <- paste0(results.folder,'gwas_effects_',trait.name,'_',suffix,'.csv')
    gwas.csv.t<- paste0(results.folder,'gwas_Tstat_',trait.name,'_',suffix,'.csv')
  }
  ############################################################

  if (reduce.K) {
    K <- reduce_kinship(K,n.pca)
  }

  ###################################################################
  # fit variance components

  # In case of unstructured (pairwise) models
  if (cov.model==2) {

    require(reshape)
    require(Matrix)

    Y.long <- melt(data=data.frame(genotype=rownames(Y), Y), variable_name = "trait.name",
                            id.vars = 'genotype',
                            measure.vars = 2:(ncol(Y)+1))

    names(Y.long)[3] <- 'pheno'


    if (run.em) {


      out <- asreml_unstructured_pairwise(d = Y.long, K=K, fix.diag=FALSE,
                                               correlation.matrix=TRUE,
                                               vE.diag=vE.diag,
                                               genotype.column=1,
                                               traitname.column=2,
                                               phenotype.column=3,
                                               covariates=integer())

      out$vG.matrix <- out$vG.matrix[colnames(Y), colnames(Y)]

      out$vE.matrix <- out$vE.matrix[colnames(Y), colnames(Y)]

      out$vG.vector <- out$vG.vector[colnames(Y)]

      out$vE.vector <- out$vE.vector[colnames(Y)]

      vG.matrix <- as.matrix(nearPD(out$vG.matrix, corr=TRUE)$mat)

      vE.matrix <- as.matrix(nearPD(out$vE.matrix, corr=TRUE)$mat)

      vG.matrix <- (matrix(sqrt(out$vG.vector)) %*% t(matrix(sqrt(out$vG.vector)))) * vG.matrix

      vE.matrix <- (matrix(sqrt(out$vE.vector)) %*% t(matrix(sqrt(out$vE.vector)))) * vE.matrix

      rownames(vG.matrix) <- rownames(vE.matrix) <- colnames(Y)

      colnames(vG.matrix) <- colnames(vE.matrix) <- colnames(Y)

      Vg <- vG.matrix

      Ve <- vE.matrix

      save(Vg, Ve,file=paste0(results.folder,varcomp.file))

    } else {

      load(file=paste0(results.folder,varcomp.file))

    }

  }

  ########################
  # In case of unstructured (pairwise) models
  if (cov.model==5) {


    if (use.heritability.package) {
      VC <- matrix(0, ncol(Y), 2)
      GBLUP <- matrix(0, nrow(Y), ncol(Y))
      for (j in 1:ncol(Y)) {
        out.h2 <- marker_h2_means(data.vector = Y[,j], geno.vector=rownames(Y), K = K, Dm=NULL, alpha = 0.05, eps = 1e-06,
                                  max.iter = 100, fix.h2 = FALSE, h2 = 0.5, grid.size=99)
        VC[j, ] <- c(out.h2$va, out.h2$ve)
        delta <- (out.h2$va / out.h2$ve)
        GBLUP[, j] <- delta * K %*% solve(delta * K + diag(nrow(Y))) %*% matrix(Y[,j])
      }
      Vg <- cov(GBLUP)
      Ve <- cov(Y - GBLUP)
      # cov2cor(Ve); cov2cor(Vg)

    } else {
      varcomp  <- EM_function_FA(Y=Y,K=K,X=data.frame(),
                            max.iter.em=max.iter.em,
                            tol.em=tol.em,
                            Cm.start=NULL,
                            Dm.start=NULL,
                            m.G=0,
                            m.E=0,
                            Cm.het=T,
                            Dm.het=T,
                            max.diag=max.diag,
                            compute.log.lik=compute.log.lik,
                            stopifdecreasing=stopifdecreasing)

      Vg <- cov(matrix(varcomp$pred$predicted,ncol=ncol(Y)))
      Ve <- cov(Y - matrix(varcomp$pred$predicted,ncol=ncol(Y)))

    }
  }
  #########################


  # In case of FA models
  if (cov.model==3) {

    if (run.em) {
      varcomp  <- EM_function_FA(Y=Y,K=K,X=X,
                          max.iter.em=max.iter.em,
                          tol.em=tol.em,
                          Cm.start=NULL,
                          Dm.start=NULL,
                          m.G=m.G,
                          m.E=m.E,
                          Cm.het=Cm.het,
                          Dm.het=Dm.het,
                          max.diag=max.diag,
                          compute.log.lik=compute.log.lik,
                          stopifdecreasing=stopifdecreasing)

      save(varcomp,file=paste0(results.folder,varcomp.file))
    } else {
      load(file=paste0(results.folder,varcomp.file))
    }

    if (snp.covariates[1]!='') {

      varcomp.red  <- EM_function_FA(Y=Y,K=K,X=X.red,
                          max.iter.em=max.iter.em,
                          tol.em=tol.em,
                          Cm.start=NULL,
                          Dm.start=NULL,
                          m.G=m.G,
                          m.E=m.E,
                          Cm.het=TRUE,
                          Dm.het=TRUE,
                          max.diag=max.diag,
                          compute.log.lik=compute.log.lik,
                          stopifdecreasing=stopifdecreasing)
    }


    Vg <- solve(varcomp$Cm)
    Ve <- solve(varcomp$Dm)

    if (snp.covariates[1]!='') {
      Vg.red <- solve(varcomp.red$Cm)
      Ve.red <- solve(varcomp.red$Dm)
    }

  }

  # In case of a user-defined Vg and Ve

  if (cov.model==4) {

    #setwd(data.folder)

    load(cov.file.name)

    stopifnot(class(Vg)=='matrix')
    stopifnot(class(Ve)=='matrix')

    stopifnot(colnames(Vg)!=NULL)
    stopifnot(colnames(Ve)!=NULL)

    stopifnot(rownames(Vg)!=NULL)
    stopifnot(rownames(Ve)!=NULL)

    stopifnot(all(colnames(Vg)==rownames(Vg)))
    stopifnot(all(colnames(Ve)==rownames(Ve)))

    stopifnot(all(colnames(Y) %in% colnames(Vg)))
    stopifnot(all(colnames(Y) %in% colnames(Ve)))

    Vg <- Vg[colnames(Y), colnames(Y)]

    Ve <- Ve[colnames(Y), colnames(Y)]

    colnames(Vg) <- rownames(Vg) <- NULL

    colnames(Ve) <- rownames(Ve) <- NULL

    #setwd(results.folder)

  }


  #######################
  # run gwas

  p  <- ncol(Y)
  n  <- nrow(X)

  nc <- ncol(X)

  if (snp.covariates[1]!='') {
    nc.red <- ncol(X.red)
  }

  w <- eigen(K)
  Dk <- diag(w$values)
  Uk <- w$vectors
  rownames.Y <- rownames(Y)

  Yt <- t(Y) %*% Uk
  colnames(Yt) <- rownames.Y

  if (nc > 0) {
    Xt <- t(X) %*% Uk
  }

  if (snp.covariates[1]!='') {
    if (nc.red > 0) {
      Xt.red <- t(X.red) %*% Uk
    }
  }

  V.inv.array <- make.V.inv.array(Vg=Vg,Ve=Ve,Dk=Dk)

  if (snp.covariates[1]!='') {
    V.inv.array.red <- make.V.inv.array(Vg=Vg.red,Ve=Ve.red,Dk=Dk)
  }

  ##########################################################################################################

  nn <- GWAS.obj$N

  means <- apply(GWAS.obj$markers[1:nn,rownames(Y)],1,mean)

  means <- means / max(GWAS.obj$markers)

  exluded.markers <- which(means < MAF | means > 1-MAF)

  if (snp.covariates[1]!='') {

    snp.covariates.numbers <- which(rownames(GWAS.obj$markers) %in% snp.covariates)
    exluded.markers <- c(exluded.markers,snp.covariates.numbers)

    extra.exluded.markers <- character()
    for (snp in snp.covariates) {
      candidates <- names(which(means==means[snp]))
      if (length(candidates) > 1) { #only the snp itself is not enough; there needs to be one other snp at least with the same maf, before proceding
        candidates.numbers <- as.integer(which(means==means[snp]))
        qwe <- apply(GWAS.obj$markers[candidates.numbers,rownames(Y)],1,function(x){identical(as.numeric(x),as.numeric(GWAS.obj$markers[snp,rownames(Y)]))})
        extra.exluded.markers <- c(extra.exluded.markers,setdiff(names(which(qwe)), snp)  )
      }
    }
    exluded.markers <- c(exluded.markers,which(rownames(GWAS.obj$markers) %in% extra.exluded.markers))
    snp.covariates.numbers <- sort(c(snp.covariates.numbers,which(rownames(GWAS.obj$markers) %in% extra.exluded.markers)))
  }

  ####################
  # scan

  if (run.gwas) {

    results <- rep(NA,nn)

    results.wald <- rep(NA,nn)

    M <- matrix(NA,nn,p)

    results[exluded.markers] <- 1

    results[exluded.markers] <- 1

    results.wald[exluded.markers] <- 1

    M[exluded.markers,] <- NA

    Tstat <- M

    if (snp.covariates[1]!='') {

      est0.red <- estimate.effects(X=Xt.red,Y=Yt,Dk=Dk,V.inv.array=V.inv.array.red,return.all.effects=T)

      fitted.mean0.red <- matrix(est0.red$effects.estimates,ncol=length(est0.red$effects.estimates)/p) %*% Xt.red

      SS0.red <- LL.quad.form.diag(Y=Yt-fitted.mean0.red,V.inv.array=V.inv.array.red)

      for (mrk in snp.covariates.numbers) {

        x <- as.numeric(t(matrix(as.numeric(GWAS.obj$markers[mrk,colnames(Yt)]))))

        xt <- t(matrix(x)) %*% Uk

        qwerty <- LRT.test(Y=Yt,X=Xt.red,x=xt,Dk=Dk,V.inv.array=V.inv.array.red,SS0=SS0.red)

        results[mrk] <- qwerty$pvalue

        results.wald[mrk] <- pchisq(sum((qwerty$effects/qwerty$effects.se)^2 ),df=p,lower.tail=F)

        M[mrk,] <- qwerty$effects

        Tstat[mrk,] <-  qwerty$effects / qwerty$effects.se

      }
    }


    est0 <- estimate.effects(X=Xt,Y=Yt,Dk=Dk,V.inv.array=V.inv.array,return.all.effects=T)

    system.time({
    fitted.mean0 <- matrix(est0$effects.estimates,ncol=length(est0$effects.estimates)/p) %*% Xt
    })

    SS0 <- LL.quad.form.diag(Y=Yt-fitted.mean0,V.inv.array=V.inv.array)

    for (mrk in setdiff(1:nn,exluded.markers)) {

      x <- as.numeric(t(matrix(as.numeric(GWAS.obj$markers[mrk,colnames(Yt)]))))

      xt <- t(matrix(x)) %*% Uk

      qwerty <- LRT.test(Y=Yt,X=Xt,x=xt,Dk=Dk,V.inv.array=V.inv.array,SS0=SS0)

      if (round(mrk/500)==(mrk/500)) {cat('Progress: ',(mrk/nn)*100,' percent','\n')}

      results[mrk] <- qwerty$pvalue

      results.wald[mrk] <- pchisq(sum((qwerty$effects/qwerty$effects.se)^2 ),df=p,lower.tail=F)

      M[mrk,] <- qwerty$effects

      Tstat[mrk,] <-  qwerty$effects / qwerty$effects.se

    }

    rownames(M) <- rownames(GWAS.obj$markers)
    colnames(M) <- colnames(Y)

    rownames(Tstat) <- rownames(GWAS.obj$markers)
    colnames(Tstat) <- colnames(Y)

    M.extended  <- data.frame(GWAS.obj$map[,c(3,1,2)],LOD_F=-log10(results),LOD_Wald=-log10(results.wald),M)
    write.csv(M.extended,quote=F,row.names=F,file=gwas.csv)

    Tstat.extended  <- data.frame(GWAS.obj$map[,c(3,1,2)],LOD_F=-log10(results),LOD_Wald=-log10(results.wald),Tstat)
    write.csv(Tstat.extended,quote=F,row.names=F,file=gwas.csv.t)

    save(M,Tstat,results,results.wald,file=paste0(gwas.file))

  } else {
    load(file=paste0(gwas.file))
  }

  ###################################################################################################
  # output :
  # results : vector of p-values
  # results.wald : vector of p-values for Wald test (for comparison; don't use it)
  # M : matrix of effect-size estimates

  #setwd(results.folder)

  library(qqman)

  mh.frame <- data.frame(BP =GWAS.obj$map$position[setdiff(1:nn,exluded.markers)],
                         CHR=GWAS.obj$map$chromosome[setdiff(1:nn,exluded.markers)],
                         P=results[setdiff(1:nn,exluded.markers)])

  ###############
  # make qq and manhattan plots, for the vector of p-values contained in results


  snps.for.qq.plot <- setdiff(1:nn,exluded.markers)


  if (snp.covariates[1]!='') {
    snps.for.qq.plot <- sort(c(snps.for.qq.plot,snp.covariates.numbers))
  }

  if (cov.model==3) {

    plot.file.name      <- paste0(results.folder,trait.name,'_',m.G,'_',m.E,suffix,'_manhattan.jpeg')
    plot.file.name.wald <- paste0(results.folder,trait.name,'_',m.G,'_',m.E,suffix,'_manhattan_wald.jpeg')

    QQ.file.name        <- paste0(results.folder,trait.name,'_',m.G,'_',m.E,suffix,'_QQplot.jpeg')
    QQ.file.name.wald   <- paste0(results.folder,trait.name,'_',m.G,'_',m.E,suffix,'_QQplot_wald.jpeg')

  } else {

    QQ.file.name        <- paste0(results.folder,trait.name,'_',suffix,'_QQplot.jpeg')
    QQ.file.name.wald   <- paste0(results.folder,trait.name,'_',suffix,'_QQplot_wald.jpeg')

    plot.file.name      <- paste0(results.folder,trait.name,'_',suffix,'_manhattan.jpeg')
    plot.file.name.wald <- paste0(results.folder,trait.name,'_',suffix,'_manhattan_wald.jpeg')

  }

  QQplotPvalues2(results[snps.for.qq.plot], main.title='QQ-plot',
                   file.name = QQ.file.name)

  QQplotPvalues2(results.wald[snps.for.qq.plot], main.title='QQ-plot',
                   file.name = QQ.file.name.wald)


  MakeLodPlotWithChromosomeColors(xvalues=GWAS.obj$map$cum.position,yvalues=-log10(results),gwas.obj=GWAS.obj,
                                           file.name=plot.file.name,jpeg.plot=T,
                                           chr.boundaries=GWAS.obj$chr.lengths.bp[-1],
                                           x.lab="Chromosomes",
                                           y.lab=expression(-log[10](p)),
                                           plot.title=trait.name,plot.type='p')

  MakeLodPlotWithChromosomeColors(xvalues=GWAS.obj$map$cum.position,yvalues=-log10(results.wald),gwas.obj=GWAS.obj,
                                           file.name=plot.file.name.wald,jpeg.plot=T,
                                           chr.boundaries=GWAS.obj$chr.lengths.bp[-1],
                                           x.lab="Chromosomes",
                                           y.lab=expression(-log[10](p)),
                                           plot.title=trait.name,plot.type='p')


}


