rm(list = ls())
library(genStatPipeline)

load("./example_data_drops/testData.Rdata")

## Create GData object
colnames(map)[1:2] <- c("chr", "pos")
pheno <- tibble::rownames_to_column(Y, var = "genotype")

gData0 <- createGData(map=map[map$chr %in% c(2,3,4,5),],
  geno=t(markers[rownames(markers) %in% rownames(map[map$chr %in% c(1,2,3,4),]),]),
  kin = K, pheno = pheno, covar = NULL)

covar <- as.data.frame(x = as.factor(sapply(strsplit(gData0$pheno[[1]]$genotype, "_"), tail, 1)))
colnames(covar) <- "region"
rownames(covar) <- gData0$pheno[[1]]$genotype
covar$altitude <- runif(nrow(covar), 0, 2000)

## For fast multitrait GWAS restrict number of traits
Y1 <- gData0$pheno[[1]][, 1:5]

## Create testdata
gDataTest <- createGData(gData0, pheno = Y1, covar = covar)
rm(list = setdiff(ls(), c("gDataTest")))

## Run single trait GWAS
testGwas1 <- runSingleTraitGwas(gData = gDataTest, thrType = 2)
## Run multi trait GWAS
testGwas2 <- runMultiTraitGwas(gData = gDataTest, subsetMarkers = TRUE, markerSubset = 40501:41500)
## Run multi trait GWAS with chromosome specific kinship matrices.
## Markersubset has to contain at least 2 chromosomes to use chr. specific kinships.
testGwas3 <- runMultiTraitGwas(gData = gDataTest,
  subsetMarkers = TRUE, markerSubset = 40501:41500, GLSMethod = 2)

## summary and plot
summary(testGwas1)
plot(testGwas1, trait = "anthesis.ARIHAS_2013_drought")
plot(testGwas1, trait = "anthesis.ARIHAS_2013_drought", type = "qq")
plot(testGwas1, type = "qtl")
