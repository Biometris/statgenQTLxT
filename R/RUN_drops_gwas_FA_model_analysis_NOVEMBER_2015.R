#
#
#

# remove all objects
rm(list=ls())
gc()

#######################################################################
# SETTINGS      (some other settings are defined when calling the functions
#                define_pheno_data and EM_function_FA below)
#######################################################################

run.em <- T    # if F, read the variance components from the file varcomp.file in results.folder

run.gwas <- T  # if F, read the gwas-results from the file gwas.file in results.folder

# the following option is still under construction; leave to zero
LOD.thr <- 0   # if larger than zero, it is assumed a GWAS was done previously with the same name
               # .. and GWAS is now only (re)run for markers with -log(p) larger than LOD.thr

MAF <- 0.05

# here are the phenotypic csv files
data.folder    <- 'D:/willem/Dropbox/drops gxe emilie/'

# here are the R-scripts
script.folder  <- 'D:/willem/research/STATISTICAL_GENETICS/AI_EM_PX_algorithms/DROPS/'

# Results will be written to this folder (needs to exist!)
results.folder <- 'D:/willem/research/STATISTICAL_GENETICS/DROPS/MULTI_Env_GWAS/'

# genotypic data (give full path; drops.RData = 50k ; otherwise drops_600k.RData)
r.image <- 'D:/willem/statistical_genetics_large_files/DROPS_data/maize_genotypic_data/drops.RData'

# select a subset of situations; to select all situations: 0
# If not 0, trait.selection should be a vector of integers, referring to situations
# in the order in which they occur in the file(s), e.g. to select the first 10, set trait.selection <- 1:10
trait.selection <- 0#1:28

# the gwas results will be stored in this RData object  (in the folder results.folder)
# if '', the name is derived from the trait name
#
#gwas.file <- 'DROPS_GWAS_GY15_23env_333k_gwas.RData'
gwas.file <- ''

# the variance components will be stored in this RData object (in the folder results.folder)
# if '', the name is derived from the trait name
#
#varcomp.file <- 'DROPS_GWAS_GY15_23env_333k_varcomp.RData'
varcomp.file <- ''

######################
# parameters for fitting the variance components
tol.em <- 0.00001

max.iter.em <- 10000

# Factor analytic models:
# Order of genetic part
m.G <- 2
# Order of environmental part
m.E <- 1

###################
# list of phenotypic files to be analyzed (one per trait)
pheno.files <- c("BLUES_Anthesis_MET.csv",
                  "BLUES_GrainNumber_MET.csv",
                  "BLUES_GrainNumberPlant_MET.csv",
                  "BLUES_GrainWeight_MET.csv",
                  "BLUES_GrainYield_MET.csv",
                  "BLUES_GrainYieldPlant_MET.csv",
                  "BLUES_PlantHeight_MET.csv",
                  "BLUES_SeedMoist_MET.csv")

#pheno.files <- c("BLUES_PlantHeight_MET.csv","BLUES_SeedMoist_MET.csv")

pheno.files <- c("grainyield_multi_env_gwas_format.csv"

########################################################################################
# START OF THE PROGRAM
########################################################################################

library(MASS)     # for the ginv function
library(mice)
library(corpcor)
library(Matrix)
library(qqman)

##############################################################################################

setwd(script.folder)

source(file="EM_functions.R")
source(file="general_functions.R")

source(file="newton_raphson_functions.R")
source(file="newton_raphson_functions_diag.R")

source(file="updateC.R")
source(file="EM_function_with_indices.R")

source(file="update_LS.R")
source(file="EM_function_LS.R")

source(file="update_FA.R")
source(file="update_FA_homogeneous_var.R")
source(file="EM_function_FA.R")

source(file='gwas_functions.R')

source(file="read_and_impute_pheno_file_FUNCTION.R")

######################################################################
# prepare data

# source(file="prepare_drops_data_emilie.R")
# source(file="read_and_impute_pheno_file.R")


for (pheno.file in pheno.files) {

  load(r.image)

  setwd(data.folder)

  pheno <- define_pheno_data(pheno.file=pheno.file,
                                standardize=T,
                                covariate.file='',
                                r.image=r.image,
                                n.impute=11,
                                trait.selection=trait.selection
                                )

  Y <- pheno$Y

  K <- pheno$K

  X <- pheno$X

  trait.name <- pheno$trait.name

  if (varcomp.file=='') {
    varcomp.file <- paste0('varcomp_',trait.name,'_',m.G,'_',m.E,'.RData')
  }

  if (gwas.file=='') {
    gwas.file <- paste0('gwas_',trait.name,'_',m.G,'_',m.E,'.RData')
  }

  ####################################################################
  # fit variance components

  if (run.em) {
    varcomp  <- EM_function_FA(Y=Y,K=K,X=X,
                        max.iter.em=max.iter.em,
                        tol.em=tol.em,
                        Cm.start=NULL,
                        Dm.start=NULL,
                        m.G=m.G,
                        m.E=m.E,
                        Cm.het=TRUE,
                        Dm.het=TRUE,
                        max.diag=100)

    save(varcomp,file=paste0(results.folder,varcomp.file))
  } else {
    load(file=paste0(results.folder,varcomp.file))
  }


  #######################
  # run gwas

  X <- matrix(rep(1,nrow(Y)))

  Vg <- solve(varcomp$Cm)
  Ve <- solve(varcomp$Dm)

  p  <- ncol(Y)
  n  <- nrow(X)
  nc <- ncol(X)

  w <- eigen(K)
  Dk <- diag(w$values)
  Uk <- w$vectors
  rownames.Y <- rownames(Y)

  Yt <- t(Y) %*% Uk
  colnames(Yt) <- rownames.Y

  if (nc > 0) {
    Xt <- t(X) %*% Uk
  }

  V.inv.array <- make.V.inv.array(Vg=Vg,Ve=Ve,Dk=Dk)

  ##########################################################################################################

  nn <- GWAS.obj$N

  means <- 0.5 * apply(GWAS.obj$markers[1:nn,rownames(Y)],1,mean)

  exluded.markers <- which(means < MAF | means > 1-MAF)


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

    est0 <- estimate.effects(X=Xt,Y=Yt,Dk=Dk,V.inv.array=V.inv.array,return.all.effects=T)

    fitted.mean0 <- matrix(est0$effects.estimates,ncol=length(est0$effects.estimates)/p) %*% Xt

    SS0 <- LL.quad.form.diag(Y=Yt-fitted.mean0,V.inv.array=V.inv.array)

    start.date <- date()

    for (mrk in setdiff(1:nn,exluded.markers)) {

      x <- as.numeric(t(matrix(as.numeric(GWAS.obj$markers[mrk,colnames(Yt)]))))

      xt <- t(matrix(x)) %*% Uk

      qwerty <- LRT.test(Y=Yt,X=Xt,x=xt,Dk=Dk,V.inv.array=V.inv.array,SS0=SS0)

      if (round(mrk/100)==(mrk/100)) {cat('Progress: ',(mrk/nn)*100,' percent','\n')}

      results[mrk] <- qwerty$pvalue

      results.wald[mrk] <- pchisq(sum((qwerty$effects/qwerty$effects.se)^2 ),df=p,lower.tail=F)

      M[mrk,] <- qwerty$effects

    }

    end.date <- date()

    save(M,results,results.wald,file=paste0(results.folder,gwas.file))

  } else {
    load(file=paste0(results.folder,gwas.file))
  }

  ###################################################################################################
  # output :
  # results : vector of p-values
  # results.wald : vector of p-values for Wald test (for comparison; don't use it)
  # M : matrix of effect-size estimates

  setwd(results.folder)

  library(qqman)

  mh.frame <- data.frame(BP =GWAS.obj$map$position[setdiff(1:nn,exluded.markers)],
                         CHR=GWAS.obj$map$chromosome[setdiff(1:nn,exluded.markers)],
                         P=results[setdiff(1:nn,exluded.markers)])

  ###############
  # make qq and manhattan plots, for the vector of p-values contained in results



  #jpeg(file=,quality=100)
  QQplotPvalues2(results[setdiff(1:nn,exluded.markers)], main.title='QQ-plot',
                 file.name = paste0(trait.name,'_QQplot.jpeg'))
  #dev.off()

  pdf(file=paste0(trait.name,'_manhattan.pdf'),width=10,height=7)
  manhattan(mh.frame)
  dev.off()

}

#########
