###
### Awesome test of the functions
###
rm(list=ls())
library(devtools)
library(tibble)

load_all(".")
## Creating object to run the analysis : createGData.R
load("./example_data_drops/testData.RData")
#Y <- rownames_to_column(Y,var="genotype")
names(map)[1:2] <- c("chr","pos")
map <- map[map$chr>0,]

# Genotyping
load("/home/mille016/Documents/PostDoc/GWAS_multitrait/Matrix_ImputedBeagle012_600K_BGA_Amaizing_Dente_Chr_All_OneMatrix_MAF05drops.Rdata")
dim(mat3)
mat4 <- mat3[,which(colnames(mat3)%in%map$snp.name)]

unique(map$chr[which(map$snp.name%in%colnames(mat3))])

# Kinship
load("/home/mille016/Documents/PostDoc/GWAS_multitrait/Kinship_Renaud_Astle_List_50K_600K_GBS_BGA_Amaizing_Dente_Chr_Panzea_.Rdata")
kin[[11]] <- NULL
kin[[1]][1:10,1:10]

Y <- read.table("~/Documents/PhD-1A/Field_Diagnostic/stress_indices/Formatting_IndexBLUES/BLUEs/BLUES_GY15.csv",sep=";",h=T)
Y <- Y[Y$Code_ID<3300,c(4,5:ncol(Y))] #
names(Y)[1] <- "genotype"
#Y <- Y[Y$genotype%in%row.names(K),]
#Y$genotype <- as.character(Y$genotype)
#Y <- droplevels(Y)

pheno <- read.csv("./example_data_drops/pheno_comparisonGWASsingletrait.csv")

#row.names(K)%in%Y$genotype
#colnames(K)%in%Y$genotype

#GData <- createGData(geno = t(markers), map = map, pheno = pheno)
GData <- createGData(geno = mat4, map = map, pheno = Y)

#rm(list = setdiff(ls(), c("GData", "K")))

## run the single trait GWAS : runSingleTraitGwas.R
GWAS <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2012_watered", environments = NULL, covar = NULL,
                          snpCovariates = NULL, K = kin, kinshipMethod = "astle", remlAlgo = 1, GLSMethod = 2,
                          sizeInclRegion = 0, useMAF = TRUE, MAF = 0.05, MAC = 10,
                          genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)

GWASkin <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2012_watered", environments = NULL, covar = NULL,
                           snpCovariates = NULL, K = NULL, kinshipMethod = "astle", remlAlgo = 1, GLSMethod = 2,
                           sizeInclRegion = 0, useMAF = TRUE, MAF = 0.05, MAC = 10,
                           genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)


GWAS2 <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2013_watered", environments = NULL, covar = NULL,
  snpCovariates = NULL, K = NULL, kinshipMethod = "GRM", remlAlgo = 1, GLSMethod = 2,
  sizeInclRegion = 0, useMAF = TRUE, MAF = 0.01, MAC = 10,
  genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)

GWAS3 <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2013_watered", environments = NULL, covar = NULL,
  snpCovariates = NULL, K = NULL, kinshipMethod = "IBS", remlAlgo = 1, GLSMethod = 2,
  sizeInclRegion = 0, useMAF = TRUE, MAF = 0.01, MAC = 10,
  genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)

GWAS4 <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2013_watered", environments = NULL, covar = NULL,
  snpCovariates = NULL, K = NULL, kinshipMethod = "vanRaden", remlAlgo = 1, GLSMethod = 2,
  sizeInclRegion = 0, useMAF = TRUE, MAF = 0.01, MAC = 10,
  genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)

GWAS5 <- runSingleTraitGwas(gData = GData, traits = "GY15.BIOGAI_2013_watered", environments = NULL, covar = NULL,
  snpCovariates = NULL, K = K, kinshipMethod = "vanRaden", remlAlgo = 1, GLSMethod = 1,
  sizeInclRegion = 0, useMAF = TRUE, MAF = 0.01, MAC = 10,
  genomicControl = FALSE, thrType = 2, alpha = 0.05, LODThr = 5, nSnpLOD = 10)

MGWAS <- runMultiTraitGwas(gData = gData, K = K, MAF = 0.01, covModel = 3)


## Comparison with older result (Sandra's pipeline, genotyping V2)
#load("~/Documents/PhD-1A/GWAS/A_50K+600K_GyGnbGsAnth/ResultsBasics/GNPasso_GY15_allARIHAS_NEW.RData")
#head(gy15)
#gy <- gy15[gy15$Situation=="GY15.BIOGAI_2013_drought"&gy15$Marker_Type=="SNP_600K",]
#gy[gy$LogPval>=5,]
load("/home/mille016/Documents/PostDoc/GWAS_multitrait/RESULTS_FASTLMM_GY15.BIOGAI_2012_drought_AllsansChrom.Rdata")
load("/home/mille016/Documents/PostDoc/GWAS_multitrait/RESULTS_FASTLMM_GY15.BIOGAI_2012_watered_AllsansChrom.Rdata")
head(resultat)


resultOrig <- read.csv("./example_data_drops/dataset_comparisonGWASsingletrait.csv")
resultOrig <- resultat
resultOrig$Marker_Name <- as.character(resultOrig$SNP)
resultNw <- GWAS$GWAResult[[1]][GWAS$GWAResult[[1]]$trait == "GY15.BIOGAI_2012_watered", ]
resultNw <- GWASkin$GWAResult[[1]][GWASkin$GWAResult[[1]]$trait == "GY15.BIOGAI_2012_watered", ]
result <- merge(resultOrig, resultNw, by.x = "Marker_Name", by.y = "snp")

jpeg(filename ="./example_data_drops/Comparison_gwasResults.jpg",width = 480, height = 480)
plot(-log10(result$Pvalue) ~ result$LOD)

abline(a = 0, b = 1, col = "red")
abline(v = 4, h = 4, col = "green",lty = 2)
dev.off()

plot(GWAS)
plot(GWAS,type="qq")


