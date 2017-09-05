context("Test single trait GWAS function")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10), matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <- paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma, pheno = list(ph1 = pheno, ph2 = pheno),
  covar = as.data.frame(covs))

result0 <- runSingleTraitGwas(gData = gDataTest, environments = 1)$GWAResult
result01 <- runSingleTraitGwas(gData = gDataTest)$GWAResult
result1 <- runSingleTraitGwas(gData = gDataTest, environments = 1, covar = "V1")$GWAResult
result2 <- runSingleTraitGwas(gData = gDataTest, environments = 1, snpCovariates = "M2")$GWAResult
result3 <- runSingleTraitGwas(gData = gDataTest, environments = 1, covar = "V1", snpCovariates = "M2")$GWAResult

test_that("runSingleTraitGwas produces output with correct dimensions", {
  expect_length(result0, 1)
  expect_length(result01, 2)
})

test_that("runSingleTraitGWas produces correct p-values", {
  expect_equal(result0[[1]]$pValue, c(0.517079439679654, 0.91343018536738, 0.628599735847542, 0.0807864803940613,
    0.857734879152352, 0.0951298087141795, 0.609379273189138, 0.999476881270353,
    0.41907977041403, 0.183886590676029, 0.973491209528092, 0.570647573548852,
    0.58845656178555, 0.367143146285207, 0.905504194974229))
  expect_equal(result1[[1]]$pValue, c(0.269946155035533, 0.985023675864513, 0.648466135580949, 0.0626938185948763,
    0.866722685001807, 0.120057200867307, 0.822217770047582, 0.896323660734102,
    0.449969869054365, 0.228248441572237, 0.946271811607473, 0.592106343634847,
    0.619273797134271, 0.405900603973107, 0.911083501168753))
  expect_equal(result2[[1]]$pValue, c(0.555483267585556, 0.946486643015959, 0.637633326342558, 0.0970517092813438,
    0.857747572789531, 0.0928128836562843, 0.668233115408406, 0.933276443599346,
    0.416819070129558, 0.212517534933692, 0.973534118714352, 0.579464414300325,
    0.648889624468317, 0.36718099462433, 0.918519084214867))
  expect_equal(result3[[1]]$pValue, c(0.30874185311605, 0.985023675864513, 0.661038060242567, 0.0830895400196558,
    0.866722685001808, 0.121327277733454, 0.832298183644918, 0.898777547575563,
    0.456921085990817, 0.264967986159878, 0.946271811607481, 0.598234089661343,
    0.657669438688024, 0.405900603973106, 0.926156990890012))
})
