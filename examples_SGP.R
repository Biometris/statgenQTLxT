rm(list = ls())
load("./example_data_drops/testData.Rdata")

colnames(map)[1:2] <- c("chr", "pos")
pheno <- tibble::rownames_to_column(Y, var = "genotype")
#gData <- createGData(map=map, geno=t(markers), kin = K, pheno = pheno, covar = NULL)

gData0 <- createGData(map=map[map$chr %in% c(2,3,4),],
  geno=t(markers[rownames(markers) %in% rownames(map[map$chr %in% c(2,3,4),]),]),
  kin = K, pheno = pheno, covar = NULL)

covar <- as.data.frame(x = as.factor(sapply(strsplit(gData0$pheno[[1]]$genotype, "_"), tail, 1)))
colnames(covar) <- "region"
rownames(covar) <- gData0$pheno[[1]]$genotype
covar$altitude <- runif(nrow(covar), 0, 2000)

# For fast multitrait GWAS
Y1 <- gData0$pheno[[1]][, 1:5]

## Create testdata
gDataTest <- createGData(gData0, pheno = Y1, covar = covar)
rm(list = setdiff(ls(), c("gDataTest")))

## Run GWAS
testGwas1 <- runSingleTraitGwas(gData = gDataTest)
testGwas2 <- runMultiTraitGwas(gData = gDataTest, subsetMarkers = TRUE, markerSubset = 40501:41500)

## summary and plot
summary(testGwas1)
plot(testGwas1, trait = "anthesis.ARIHAS_2013_drought")
plot(testGwas1, trait = "anthesis.ARIHAS_2013_drought", type = "qq")

