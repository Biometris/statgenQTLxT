
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
