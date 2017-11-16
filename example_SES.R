rm(list = ls())
load("D:/R packages/SESvdH/genomic_selection/2016/GPdata_imputed_all.RData")
map <- markermap[, c("chr", "pos")]
rownames(map) <- markermap$name
markers <- data[which(colnames(data) == "AX124322256"):ncol(data)]
rownames(markers) <- data$Seedname
markers <- markers[data$Phenotyped == "y",]
pheno <- data[, which(colnames(data) == "THA"):which(colnames(data) == "mMN")]
pheno <- data.frame(genotype = data$Seedname, pheno, stringsAsFactors = FALSE)
covar <- data[, c("Pool", "Fam", "PepiNo", "Phenotyped")]
rownames(covar) <- data$Seedname
gDataSes <- createGData(map = map, geno = markers, pheno = pheno, covar = covar)
rm(list = setdiff(ls(), "gDataSes"))

GWASes <- runSingleTraitGwas(gData = gDataSes,
  traits = NULL,
  environments =  NULL,
  covar = NULL,
  snpCovariates = NULL,
  K = NULL,
  remlAlgo = 1,
  GLSMethod = 2,
  sizeInclRegion = 0,
  minR2 = 0.5,
  useMAF = TRUE,
  MAF = 0,
  MAC = 10,
  genomicControl = FALSE,
  thrType = 1,
  alpha = 0.005,
  LODThr = 4)

GWASorig <- read.csv("D:/R packages/SESvdH/genomic_selection/2016/GWAS_analyses.csv")
plot(x = GWASorig$minlog10p, y = -log10(GWASes$GWAResult[[1]]$pValue),
  col = adjustcolor(1:9, alpha.f = 0.5))
legend(40, 25, unique(GWASorig$chr), col = 1:length(GWASorig$chr), pch = 1)
abline(0,1, col = "blue")

plot(x = GWASorig$minlog10p, y = -log10(GWASes$GWAResult[[1]]$pValue),
  col = adjustcolor(1:8, alpha.f = 0.3))
legend(33, 25, unique(GWASorig$trait), col = 1:length(GWASorig$trait), pch = 1)
abline(0,1, col = "blue")

