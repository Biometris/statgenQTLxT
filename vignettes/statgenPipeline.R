## ----setup, include = FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)
library(statgenPipeline)
options(width = 100, digits = 2)

## ----loadData-------------------------------------------------------------------------------------
library(statgenGWAS)
data(dropsMarkers)
data(dropsMap)
data(dropsPheno)

## ----convertMarkers-------------------------------------------------------------------------------
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## ----convertMap-----------------------------------------------------------------------------------
## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

## ----createGdata----------------------------------------------------------------------------------
## Create a gData object containing map and marker information.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap)

## ----addPheno-------------------------------------------------------------------------------------
## Convert phenotypic data to a list.
dropsPhenoList <- split(x = dropsPheno, f = dropsPheno[["Experiment"]])
## Rename Variety_ID to genotype and select relevant columns.
dropsPhenoList <- lapply(X = dropsPhenoList, FUN = function(trial) {
  colnames(trial)[colnames(trial) == "Variety_ID"] <- "genotype"
  trial <- trial[c("genotype", "grain.yield", "grain.number", "anthesis", 
                   "silking", "plant.height", "ear.height")]
  return(trial)
})
## Add phenotypic data to gDataDrops.
gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList)

## ----sumGData-------------------------------------------------------------------------------------
## Summarize gDataDrops.
summary(gDataDrops, trials = "Mur13W")

## ----removeDupMarkers-----------------------------------------------------------------------------
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 

## ----addMissings----------------------------------------------------------------------------------
## Copy gData object.
gDataDropsMiss <- gDataDrops
## Add random missing values to 1% of the values in the marker matrix.
set.seed(1)
nVal <- nrow(gDataDropsMiss$markers) * ncol(gDataDropsMiss$markers)
gDataDropsMiss$markers[sample(x = 1:nVal, size = nVal / 100)] <- NA

## ----imputeMissings-------------------------------------------------------------------------------
## Impute missing values with random value.
## Remove SNPs and genotypes with proportion of NA larger than 0.01.
gDataDropsImputed <- codeMarkers(gData = gDataDropsMiss, 
                                 nMissGeno = 0.01, 
                                 nMiss = 0.01, 
                                 impute = TRUE, 
                                 imputeType = "random", 
                                 verbose = TRUE)

## ----imputeMissingsBeagle, eval=FALSE-------------------------------------------------------------
#  ## Impute missing values using beagle software.
#  gDataDropsImputedBeagle <- codeMarkers(gData = gDataDropsMiss,
#                                         impute = TRUE,
#                                         imputeType = "beagle",
#                                         verbose = TRUE)

## ----mtg------------------------------------------------------------------------------------------
## Run multi trait GWAS for traits 'grain.yield' and 'anthesis' for trial Mur13W.
GWASDrops <- runMultiTraitGwas(gData = gDataDropsDedup, trials = "Mur13W", covModel = "fa")

## ----gwaRes---------------------------------------------------------------------------------------
head(GWASDrops$GWAResult$Mur13W)

## ----signSnp--------------------------------------------------------------------------------------
print(GWASDrops$signSnp$Mur13W, row.names = FALSE)

## ----sumMtg---------------------------------------------------------------------------------------
## Create summary of GWASDrops.
summary(GWASDrops)

## ----qqMtg----------------------------------------------------------------------------------------
## Plot a qq plot of GWAS Drops.
plot(GWASDrops, plotType = "qq")

## ----manhattanMtg---------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
plot(GWASDrops, plotType = "manhattan")

## ----manhattanMtgThr------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
## Set significance threshold to 4 and only plot chromosomes 6 to 8.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", yThr = 4, chr = 6:8)

## ----manhattanLod---------------------------------------------------------------------------------
## Plot a manhattan plot of GWAS Drops.
## Plot only 5% of SNPs with a LOD below 3.
set.seed(1)
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", lod = 3)

## ----qtlMtgNorm-----------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 3 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 3, normalize = TRUE)

## ----mtgChrSpec-----------------------------------------------------------------------------------
## Run multi trait GWAS for trial 'Mur13W'.
## Use chromosome specific kinship matrices computed using method of van Raden.
GWASDropsChrSpec <- runMultiTraitGwas(gData = gDataDropsDedup, 
                                      trials = "Mur13W",
                                      GLSMethod = "multi",
                                      kinshipMethod = "vanRaden",
                                      covModel = "fa")

## ----mtgSNPFixThr---------------------------------------------------------------------------------
## Run multi trait GWAS for Mur13W.
## Use a fixed significance threshold of 4.
GWASDropsFixThr <- runMultiTraitGwas(gData = gDataDropsDedup,
                                     trials = "Mur13W", 
                                     covModel = "fa")

## ----mtgSNPNR-------------------------------------------------------------------------------------
## Run multi trait GWAS for for Mur13W.
## Use a factor analytic model for computing the variance components.
GWASDropsFA <- runMultiTraitGwas(gData = gDataDropsDedup,
                                 trials = "Mur13W",
                                 covModel = "fa")

## Rerun the analysis, using the variance components computed in the 
## previous model as inputs.
GWASDropsFA2 <- runMultiTraitGwas(gData = gDataDropsDedup,
                                  trials = "Mur13W",
                                  fitVarComp  = FALSE,
                                  Vg = GWASDropsFA$GWASInfo$varComp$Vg,
                                  Ve = GWASDropsFA$GWASInfo$varComp$Ve)

## ----mtgSNPCovar----------------------------------------------------------------------------------
## Run multi trait GWAS for Mur13W.
## Use PZE-106021410, the most significant SNP, a SNP covariate.
GWASDropsSnpCov <- runMultiTraitGwas(gData = gDataDropsDedup,
                                     trials = "Mur13W",
                                     snpCov = "PZE-106021410",
                                     covModel = "fa")

## ----mtgMAF---------------------------------------------------------------------------------------
## Run multi trait GWAS for Mur13W.
## Only include SNPs that have a MAF of 0.05 or higher.
GWASDropsMAF <- runMultiTraitGwas(gData = gDataDropsDedup,
                                  trials = "Mur13W",
                                  covModel = "fa",
                                  MAF = 0.05)

## ----mtgCommon------------------------------------------------------------------------------------
## Run multi trait GWAS for Mur13W.
## Fit an additional common sNP effect model.
GWASDropsCommon <- runMultiTraitGwas(gData = gDataDropsDedup,
                                     trials = "Mur13W",
                                     covModel = "fa",
                                     estCom = TRUE)
head(GWASDropsCommon$GWAResult$Mur13W)


