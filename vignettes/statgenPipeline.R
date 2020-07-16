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

## ----convertMarkers, echo=FALSE-------------------------------------------------------------------
## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## ----convertMap, echo=FALSE-----------------------------------------------------------------------
## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

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

## ----createGdata----------------------------------------------------------------------------------
## Create a gData object containing map, marker and phenotypic information.
gDataDrops <- createGData(geno = dropsMarkers,
                          map = dropsMap, 
                          pheno = dropsPhenoList)

## ----sumGData-------------------------------------------------------------------------------------
## Summarize gDataDrops.
summary(gDataDrops, trials = "Mur13W")

## ----removeDupMarkers-----------------------------------------------------------------------------
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 

## ----mtg------------------------------------------------------------------------------------------
## Run multi-trait GWAS for xx traits in trial Mur13W.
GWASDrops <- runMultiTraitGwas(gData = gDataDropsDedup, 
                               trials = "Mur13W", 
                               covModel = "fa")

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

## ----qtlMtgNorm-----------------------------------------------------------------------------------
## Plot a qtl plot of GWAS Drops for Mur13W.
## Set significance threshold to 4 and normalize effect estimates.
plot(GWASDrops, plotType = "qtl", yThr = 4, normalize = TRUE)

## ----mtgChrSpec, eval = FALSE---------------------------------------------------------------------
#  ## Run multi-trait GWAS for trial 'Mur13W'.
#  ## Use chromosome specific kinship matrices computed using method of van Raden.
#  GWASDropsChrSpec <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                        trials = "Mur13W",
#                                        GLSMethod = "multi",
#                                        kinshipMethod = "vanRaden",
#                                        covModel = "fa")

## ----addPhenoxE-----------------------------------------------------------------------------------
## Reshape phenotypic data to data.frame in wide format containing only grain yield.
PhenoDat <- reshape(dropsPheno[,c("Experiment","Variety_ID","grain.yield")], 
                    timevar = "Experiment", 
                    idvar = "Variety_ID", 
                    direction = "wide", 
                    v.names = "grain.yield")
## Rename Variety_ID to genotype and other columns with the experiment name only.
names(PhenoDat)[1] <- "genotype"
names(PhenoDat)[2:ncol(PhenoDat)] <-  gsub("grain.yield.","",names(PhenoDat)[2:ncol(PhenoDat)] )

## ----createGdataxE--------------------------------------------------------------------------------
## Create a gData object containing map, marker and phenotypic information.
gDataDropsxE <- createGData(geno = dropsMarkers,
                            map = dropsMap, 
                            pheno = PhenoDat)
summary(gDataDropsxE)

## ----removeDupMarkersxE---------------------------------------------------------------------------
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedupxE <- codeMarkers(gDataDropsxE, impute = FALSE, verbose = TRUE) 

## ----mtgxE----------------------------------------------------------------------------------------
## Run multi-trait GWAS for xx traits in trial Mur13W.
GWASDropsxE <- runMultiTraitGwas(gData = gDataDropsDedupxE, 
                                 covModel = "fa")

## ----gwaResxE-------------------------------------------------------------------------------------
head(GWASDropsxE$GWAResult$PhenoDat)

## ----signSnpxE------------------------------------------------------------------------------------
head(GWASDropsxE$signSnp$PhenoDat, row.names = FALSE)

## ----sumMtgxE-------------------------------------------------------------------------------------
summary(GWASDropsxE)

## ----qqMtgxE--------------------------------------------------------------------------------------
plot(GWASDropsxE, plotType = "qq")

## ----manhattanMtgxE-------------------------------------------------------------------------------
plot(GWASDropsxE, plotType = "manhattan")

## ----qtlMtgNormxE---------------------------------------------------------------------------------
## Set significance threshold to 6 and do not normalize effect estimates.
plot(GWASDropsxE, plotType = "qtl", yThr = 6, normalize = FALSE)

## ----mtgSNPFixThr, eval = FALSE-------------------------------------------------------------------
#  ## Run multi-trait GWAS for Mur13W.
#  ## Use a fixed significance threshold of 4.
#  GWASDropsFixThr <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                       trials = "Mur13W",
#                                       covModel = "fa")

## ----mtgSNPNR, eval = FALSE-----------------------------------------------------------------------
#  ## Run multi-trait GWAS for for Mur13W.
#  ## Use a factor analytic model for computing the variance components.
#  GWASDropsFA <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                   trials = "Mur13W",
#                                   covModel = "fa")
#  
#  ## Rerun the analysis, using the variance components computed in the
#  ## previous model as inputs.
#  GWASDropsFA2 <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                    trials = "Mur13W",
#                                    fitVarComp  = FALSE,
#                                    Vg = GWASDropsFA$GWASInfo$varComp$Vg,
#                                    Ve = GWASDropsFA$GWASInfo$varComp$Ve)

## ----mtgSNPCovar, eval = FALSE--------------------------------------------------------------------
#  ## Run multi-trait GWAS for Mur13W.
#  ## Use PZE-106021410, the most significant SNP, a SNP covariate.
#  GWASDropsSnpCov <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                       trials = "Mur13W",
#                                       snpCov = "PZE-106021410",
#                                       covModel = "fa")

## ----mtgMAF, eval = FALSE-------------------------------------------------------------------------
#  ## Run multi-trait GWAS for Mur13W.
#  ## Only include SNPs that have a MAF of 0.05 or higher.
#  GWASDropsMAF <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                    trials = "Mur13W",
#                                    covModel = "fa",
#                                    MAF = 0.05)

## ----mtgCommon, eval = FALSE----------------------------------------------------------------------
#  ## Run multi-trait GWAS for Mur13W.
#  ## Fit an additional common sNP effect model.
#  GWASDropsCommon <- runMultiTraitGwas(gData = gDataDropsDedup,
#                                       trials = "Mur13W",
#                                       covModel = "fa",
#                                       estCom = TRUE)
#  head(GWASDropsCommon$GWAResult$Mur13W)
#  

