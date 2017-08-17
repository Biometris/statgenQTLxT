###
### Awesome test of the functions
###
rm(list=ls())
library(devtools)
library(tibble)

load_all()
## Creating object to run the analysis : createGData.R
#load("./example_data_drops/drops_600k.RData")
load("./example_data_drops/testData.RData")
#Y <- rownames_to_column(Y,var="genotype")
names(map)[1:2] <- c("chr","pos")

Y <- read.table("~/Documents/PhD-1A/Field_Diagnostic/stress_indices/Formatting_IndexBLUES/BLUEs/BLUES_GY15.csv",sep=";",h=T)
Y <- Y[Y$Code_ID<3300,c(3,5:ncol(Y))] #
names(Y)[1] <- "genotype"
Y <- Y[Y$genotype%in%row.names(K),]
Y$genotype <- as.character(Y$genotype)
Y <- droplevels(Y)

#row.names(K)%in%Y$genotype
#colnames(K)%in%Y$genotype

GData <- createGData(geno=t(markers), map=map, kin=K, pheno=Y)

## run the single trait GWAS : runSingleTraitGwas.R
GWAS <- runSingleTraitGwas(gData=GData,traits = 4, fields = NULL, covar = NULL,
                          snpCovariates = NULL, K = NULL, remlAlgo = 1, GLSMethod = 2,
                          sizeInclRegion = 0, minR2, useMAF = TRUE, MAF = 0.05, MAC = 10,
                          genomicControl = FALSE, thrType = 2, alpha = 0.05 , LODThr = 5 ,nSnpLOD = 10)


## Comparison with older result (Sandra's pipeline, genotyping V2)
load("~/Documents/PhD-1A/GWAS/A_50K+600K_GyGnbGsAnth/ResultsBasics/GNPasso_GY15_allARIHAS_NEW.RData")
head(gy15)
gy <- gy15[gy15$Situation=="GY15.BIOGAI_2013_drought"&gy15$Marker_Type=="SNP_600K",]
gy[gy$LogPval>=5,]

bartjan <- GWAS$GWAResult[[1]]

result <- merge(gy,bartjan,by.x="Marker_Name",by.y="snp")
jpeg(filename ="~/statGenPipeline/example_data_drops/Comparison_gwasResults.jpg",width = 480, height = 480)
plot(result$LogPval~result$LOD)
abline(v=3,col='red',lty=2)
abline(h=3,col='red',lty=2)
abline(v=4,col='green',lty=2)
abline(h=4,col='green',lty=2)
dev.off()



