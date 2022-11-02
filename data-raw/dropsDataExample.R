## Create a gData object based on full drops data with 3 traits and 10% of the
## markers.

library(statgenGWAS)

## Add genotypes as row names of dropsMarkers and drop Ind column.
set.seed(1234)
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[-1]
dropsMarkers <- dropsMarkers[, sample(x = seq_len(ncol(dropsMarkers)),
                                      size = ncol(dropsMarkers) / 10)]

## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
colnames(dropsMap)[2:3] <- c("chr", "pos")
## Restrict to selected markers.
dropsMap <- dropsMap[rownames(dropsMap) %in% colnames(dropsMarkers), ]

## Convert phenotypic data to a list.
colnames(dropsPheno)[1] <- "genotype"
dropsPheno <- dropsPheno[c("Experiment", "genotype", "grain.yield",
                           "grain.number", "anthesis")]
dropsPhenoList <- split(x = dropsPheno, f = dropsPheno[["Experiment"]])
dropsPhenoList <- lapply(dropsPhenoList, FUN = `[`, -1)

## Create a gData object containing map, marker and phenotypic information.
gDataDropsRestr <- createGData(geno = dropsMarkers,
                               map = dropsMap,
                               pheno = dropsPhenoList)

## Remove duplicate SNPs from gDataDrops.
gDataDropsRestr <- codeMarkers(gDataDropsRestr, impute = FALSE, verbose = TRUE)

usethis::use_data(gDataDropsRestr, overwrite = TRUE)

