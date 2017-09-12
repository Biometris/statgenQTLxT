load('./example_data_drops/atwell_regmap.RData')

map <- GWAS.obj$map[, 1:2]
colnames(map) <- c('chr', 'pos')
rownames(map) <- rownames(GWAS.obj$markers)

pheno <- matrix(rnorm(n = 3 * GWAS.obj$n), ncol = 3)
pheno <- data.frame(genotype = GWAS.obj$plant.names, pheno)
names(pheno)[-1] <- paste0('trait', 1:3)

g <- createGData(geno = t(GWAS.obj$markers), map = map, pheno = pheno)

rm(list = setdiff(ls(), g))


test.gwas1 <- runSingleTraitGwas(gData = g)
# function unknown !

test.gwas2 <- runMultiTraitGwas(gData = g)
#  Error in markersRed[rownames(Y), rownames(mapRed)] :
# subscript out of bounds


