context("STG IBD Helper functions")

set.seed(1234)
y <- 1:10
X <- matrix(runif(n = 90), nrow = 30)
X <- X / rowSums(X)
dim(X) <- c(10, 3, 3)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- paste0("M", 1:3)
dimnames(X)[[3]] <- paste0("A", 1:3)
gDataTest <- createGData(geno = X, kin = Sigma, pheno = pheno,
                         covar = as.data.frame(covs))

GLS0 <- fastGLSIBD(y = y, X = X, Sigma = Sigma, ref = 1)
test_that("fastGLSIBD produces correct output structure", {
  expect_is(GLS0, "matrix")
  expect_equal(dim(GLS0), c(3, 4))
  expect_equal(rownames(GLS0), paste0("M", 1:3))
  expect_equal(colnames(GLS0), c("pValue", "RLR2", "A2", "A3"))
})

test_that("fastGLSIBD produces correct output value", {
  expect_equivalent(GLS0,
                    c(0.151983822698907, 0.0105999550089082, 0.152511992087561,
                      0.999999999999998, 1, 0.999999999999998, 13.7799500899171,
                      24.7422182112529, -5.20527153225294, 8.87407632732625,
                      5.01580069441162, -14.6762391396771))
})

test_that("choice of ref allele doesn't affect fastGLSIBD output", {
  GLS1 <- fastGLSIBD(y = y, X = X, Sigma = Sigma, ref = 2)
  expect_equal(GLS0[, "pValue"], GLS1[, "pValue"])
  expect_equal(GLS0[, "A2"], -GLS1[, "A1"])
})

test_that("fastGLSIBD with covariates produces correct output value", {
  GLS2 <- fastGLSIBD(y = y, X = X, Sigma = Sigma, covs = covs, ref = 1)
  expect_equivalent(GLS2,
                    c(0.537681841019699, 0.430915845432764, 0.188623090467351,
                      0.982210259336496, 0.994705294971946, 0.99986695817212,
                      7.28929606576913, 13.6625002952502, -1.1456663141888,
                      -3.41870440778646, 9.27826137792016, -15.4883014514782))
})

