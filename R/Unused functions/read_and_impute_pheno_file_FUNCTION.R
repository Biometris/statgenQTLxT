##########################################################
# # Phenotypic data
#
# # Vector of phenotypic files to be analyzed (e.g. one per trait)
# #pheno.files <- c('grainnumber_multi_env_gwas_dec2015.csv','anthesis_multi_env_gwas_dec2015.csv')
# #
# # Each file needs to be a csv file (comma-separated; not with semi-colons !),
# # and needs to contain genotypic means for a number of environments or traits
# # The first column(s) (there may be more than 1) should contain genotype labels.
# # In case of multi-trait data, each file should contain the set of traits that should
# # be analyzed together.
# # In case of multi-environment data, each file should contain all environments for one trait.
# # The scripts are not really made for multi-trait/environment simultaneously, although...
# pheno.files <- './example_data_DROPS/anthesis_multi_env_gwas_dec2015.csv'
#
# # In the phenotypic files, do all columns start with the trait names ? (as in DROPS)
# read.trait.names.from.file <- TRUE
#
# # If read.trait.names.from.file==F :
# # Should correspond to pheno.files
# # Is only used to create the names of the output files
# # Can be a vector of '', e.g. c('','','') in case pheno.files has length 3
# # Only important in case you do a multi-environment GWAS for several traits separately,
# #  pheno.files containing the files for each trait
# trait.names <- ''
#
# # The first n.omit columns: these are all columns not containing phenotypic data,
# # in particular genotype labels
# n.omit <- 1
#
# # The column-name containing the genotype/accession/hybrid etc labels.
# #   Should be one of the first n.omit (see above) columns
# #   Should match with the names in GWAS.obj$plant.names , contained in r.image (see above)
# #   Always specify; cannot be omitted if n.omit == 1
# geno.ID <- 'Accession_ID'
#
# # select a subset of situations; to select all situations: 0
# # If not 0, trait.selection should be a vector of integers, referring to situations
# # in the order in which they occur in the file(s), e.g. to select the first 10, set trait.selection <- 1:10
# trait.selection <- 0
#
# # Standardize all phenotypic columns ?
# standardize <- F
#
# # Should imputation be performed ? Can only be F in case of balanced trials
# no.imputation <- TRUE
# for (pheno.file in pheno.files) {
#   #pheno.file = pheno.files[1]



define_pheno_data <- function(pheno.file,
  standardize=T,
  covariate.file='',
  snp.covariates='',
  r.image='drops.RData',
  n.impute=11,
  no.imputation=FALSE,
  trait.selection=integer(),
  n.omit=4,
  geno.ID='Accession_ID',
  read.trait.name.from.file=TRUE,
  trait.name=NULL
) {

  # pheno.file : should not contain duplicate hybrids, e.g. B73

  # standardize : standardize each environment ?

  # covariate.file : if no covariates are to be included, put covariate.file <- ''

  # n.impute :
  #   impute missing values for hybrids for which the number of
  #   missing situations is at most n.impute; otherwise omit these hybrids

  # r.image : R-image of the genotypic data


  # n.impute=11; no.imputation=T

  ################################################################

  if (!no.imputation) {
    miceImpute <- function(X) {
      require(mice)
      X.imp <- mice(X)
      X.imp <- complete(X.imp)
      return(X.imp)
    }
  }

  a <- read.csv(file=pheno.file)

  b <- a[,-(1:n.omit)]

  #rownames(b) <- a$Accession_ID
  rownames(b) <- a[,geno.ID]

  if (sum(trait.selection==0)) {trait.selection <- 1:ncol(b)}

  b <- b[apply(b,1,function(x){sum(is.na(x))}) <= n.impute,]

  load(r.image)

  b <- b[which(rownames(b) %in% GWAS.obj$plant.names),]


  # impute
  if (!no.imputation) {
    b <- miceImpute(b)
  }

  # subset of  traits

  b.rnames <- rownames(b)
  b.cnames <- colnames(b)
  b <- as.data.frame(b[,trait.selection])
  rownames(b) <- b.rnames
  colnames(b) <- b.cnames[trait.selection]

  Y <- b

  Y.unstandardized <- Y

  if (standardize) {
    for (j in 1:ncol(Y)) {Y[,j] <- scale(Y[,j])}
  }

  K <- GWAS.obj$kinship[rownames(b), rownames(b)]

  #  if (covariate.file=='') {
  #    X <- matrix(1,nrow(Y),1)
  #    rownames(X) <- rownames(Y)
  #  } else {
  #    # TO DO
  #  }

  if (covariate.file=='') {
    if (standardize) {
      X <- matrix(1,nrow(Y),0)
    } else {
      X <- matrix(1,nrow(Y),1)
    }
    rownames(X) <- rownames(Y)
  } else {
    # TO DO
  }

  if (snp.covariates[1]!='') {
    stopifnot(all(snp.covariates %in% rownames(GWAS.obj$markers)))
    #for (snp in snp.covariates) {
    #}
    X <- cbind(X,t(as.matrix(GWAS.obj$markers[snp.covariates,rownames(X)])))
  }

  if (read.trait.name.from.file==TRUE) {
    trait.name <- unlist(strsplit(colnames(Y), "[.]"))[1]
  } else {
    # In this case, trait.name should have been given as input
    stopifnot(!is.null(trait.name))
  }



  # output : Y,K,X
return(list(Y=Y,K=K,X=X,trait.name=trait.name,Y.unstandardized=Y.unstandardized))
}
