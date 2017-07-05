
# standardize each environment ?
standardize <- F

# should not contain duplicate hybrids, e.g. B73
pheno.file <- 'BLUES_Anthesis_MET.csv'
# D:/willem/research/STATISTICAL_GENETICS/DROPS/files_emilie_19june_2015/

covariate.file <- '' # if no covariates are to be included, put covariate.file <- ''

# 50 k
r.image   <- 'D:/willem/statistical_genetics_large_files/DROPS_data/maize_genotypic_data/drops.RData'
# 333 k
#r.image   <- 'D:/willem/statistical_genetics_large_files/DROPS_data/maize_genotypic_data/drops_600k.RData'

n.impute <- 11 # impute missing values for hybrids for which the number of
               # missing situations is at most n.impute; otherwise omit these hybrids

trait.selection <- 1:25

#50k    <- TRUE # if F, use the 333k set

################################################################

miceImpute <- function(X) {
  require(mice)
  X.imp <- mice(X)
  X.imp <- complete(X.imp)
return(X.imp)
}


a <- read.csv(file=pheno.file)

b <- a[,-(1:4)]
rownames(b) <- a$Accession_ID

b <- b[apply(b,1,function(x){sum(is.na(x))}) <= n.impute,]

if (standardize) {
  for (j in 1:ncol(b)) {b[,j] <- scale(b[,j])}
}

load(r.image)

b <- b[which(rownames(b) %in% GWAS.obj$plant.names),]

# impute
b <- miceImpute(b)

# subset of  traits
b <- b[,trait.selection]

Y <- b

K <- GWAS.obj$kinship[rownames(b), rownames(b)]

if (covariate.file=='') {
  X <- matrix(1,nrow(Y),1)
  rownames(X) <- rownames(Y)
}

# output : Y,K,X











