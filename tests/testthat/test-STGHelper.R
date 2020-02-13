context("STG Helper functions")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma, pheno = pheno,
                         covar = as.data.frame(covs))

test_that("function exclMarkers functions properly", {
  markers <- matrix(c(0, 1, 0, 1, 2, 1, 0, 1, 0, 2, 1, 2), ncol = 4,
                    dimnames = list(paste0("IND", 1:3), paste0("SNP", 1:4)))
  allFreq <- colMeans(markers, na.rm = TRUE)
  expect_length(exclMarkers(snpCov = NULL, markers = markers,
                            allFreq = allFreq), 0)
  expect_equal(exclMarkers(snpCov = "SNP2", markers = markers,
                           allFreq = allFreq), 2)
  expect_equal(exclMarkers(snpCov = "SNP1", markers = markers,
                           allFreq = allFreq), c(1, 3))
  expect_equal(exclMarkers(snpCov = c("SNP1", "SNP3"), markers = markers,
                           allFreq = allFreq), c(1, 3))
  expect_equal(exclMarkers(snpCov = c("SNP1", "SNP2", "SNP3"),
                           markers = markers, allFreq = allFreq), c(1, 3, 2))
})

test_that("genCtrlPVals produces correct new p-values", {
  expect_equal(genCtrlPVals(pVals = .5, nObs = 10)[[1]], 0.5)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[1]],
               c(0.410638105484779, 0.63432328826532))
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[1]],
               c(0.410588927034021, 0.629472060364479))
})

test_that("genCtrlPVals produces correct inflation factor", {
  expect_equal(genCtrlPVals(pVals = .5, nObs = 10)[[2]], 1)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 10)[[2]],
               2.04152779518634)
  expect_equal(genCtrlPVals(pVals = c(0.25, 0.5), nObs = 1e6)[[2]],
               1.95438310682772)
})

