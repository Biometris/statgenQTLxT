context("Test helper functions for single trait GWAS.")

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
gDataTest <- createGData(map = map, geno = X, kin = Sigma, pheno = pheno, covar = as.data.frame(covs))

test_that("runEmma produces correct results with default settings", {
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1)[[1]]),
    c(0.0214570264923801, 1.85343825287729))
  expect_equal(unname(runEmma(gData = gDataTest, trait = "X1", environment = 1)[[1]]),
    c(0.0214570264923801, 1.85343825287729))
  expect_equal(unname(runEmma(gData = gDataTest, trait = "X1", environment = "pheno")[[1]]),
    c(0.0214570264923801, 1.85343825287729))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1, K = Sigma)[[1]]),
    c(0.0214570264923801, 1.85343825287729))
})

test_that("runEmma produces correct results with covariates", {
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1, covar = 1)[[1]]),
    c(8.76962021955844e-05, 1.93163739799549))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1, covar = "V1")[[1]]),
    c(8.76962021955844e-05, 1.93163739799549))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1, K = Sigma, covar = 1)[[1]]),
    c(8.76962021955844e-05, 1.93163739799549))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1, snpName = "M1")[[1]]),
    c(0.199457967788946, 1.82420297857926))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1,
    K = Sigma, snpName = "M1")[[1]]), c(0.199457967788946, 1.82420297857926))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1,
    covar = 1, snpName = "M1")[[1]]), c(0.140158216807375, 1.71023973698722))
  expect_equal(unname(runEmma(gData = gDataTest, trait = 2, environment = 1,
    K = Sigma, covar = 1, snpName = "M1")[[1]]), c(0.140158216807375, 1.71023973698722))
})

test_that("extra options in runEmma don't significantly change results", {
  ## Compute base value
  runEmma0 <- runEmma(gData = gDataTest, trait = 2, environment = 1)[[1]]
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, nGrids = 50)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, nGrids = 500)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, lLim = -100)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, lLim = -20)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, uLim = 20)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, uLim = 100)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, eps = 1e-5)[[1]],
    runEmma0, tolerance = 1e-6)
  expect_equal(runEmma(gData = gDataTest, trait = 2, environment = 1, nGrids = 500,
    lLim = -100, uLim = 100)[[1]], runEmma0, tolerance = 1e-6)
})

test_that("emmaEigenR produces correct results", {
  expect_equal(dim(emmaEigenR(K = Sigma, X = covs)[[2]]), c(10, 8))
  expect_equal(emmaEigenR(K = Sigma[1:3, 1:3], X = covs[1:3, ])[[1]], 1.96517890694063)
  expect_equal(emmaEigenR(K = Sigma[1:3, 1:3], X = covs[1:3, ])[[2]],
    c(0.118208501113814, -0.364967722967346, 0.923485414858543))
  expect_error(emmaEigenR(K = Sigma, X = covs[1:3, ]))
  expect_error(emmaEigenR(K = Sigma[1:3, ], X = covs))
})

test_that("EMMAREMLLL produces correct output", {
  expect_equal(emmaREMLLL(logDelta = -1, lambda = 1, etas1 = 2, n = 0, t = 0, etas2 = 0), -2.11208571376462)
  expect_equal(emmaREMLLL(logDelta = -1, lambda = 1:10, etas1 = 1:10, n = 0, t = 0, etas2 = 0), -30.4461091155684)
  expect_equal(emmaREMLLL(logDelta = -1, lambda = 1:10, etas1 = 1:10, n = 3, t = 2, etas2 = 1:2),
    c(-31.9439236241564, -32.2122231996107))
})

test_that("fastGLS without covariates procudes correct output", {
  GLS0 <- fastGLS(y = y, X = X, Sigma = Sigma)
  expect_equal(GLS0[, 1], c(0.191990244479038, 0.0346367487131218, 0.297099155797429))
  expect_equal(GLS0[, 2], c(-0.765513606856416, -2.9009783711453, 0.760762479686889))
  expect_equal(GLS0[, 3], c(0.125433669816756, 0.319994381086895, 0.152892637284377))
  expect_equal(GLS0[, 4], c(0.975876824782575, 0.999730440588029, 0.91590885222019))
})

test_that("fastGLS with covariates procudes correct output", {
  GLS1 <- fastGLS(y = y, X = X, Sigma = Sigma, covs = covs)
  expect_equal(GLS1[, 1], c(0.729670715779269, 0.0632229836715346, 0.489762145590089))
  expect_equal(GLS1[, 2], c(0.779388868965725, -3.50930856822299, 1.51419592119873))
  expect_equal(GLS1[, 3], c(0.483570573233199, 0.467906979316023, 0.47775066615577))
  expect_equal(GLS1[, 4], c(0.228770901452724, 0.996393508816207, 0.633782123371107))
})

test_that("fastGLS is independent of dimensions and number of chunks", {
  expect_equal(fastGLS(y = y, X = X, Sigma = Sigma, nChunks = 5),
    fastGLS(y = y, X = X, Sigma = Sigma, nChunks = 2))
  expect_equal(fastGLS(y = y, X = X[, 1:2], Sigma = Sigma, nChunks = 5),
    fastGLS(y = y, X = X, Sigma = Sigma, nChunks = 2)[1:2, ])
  expect_equal(fastGLS(y = y, X = X, Sigma = Sigma, covs = covs, nChunks = 5),
    fastGLS(y = y, X = X, Sigma = Sigma, covs = covs, nChunks = 2))
  expect_equal(fastGLS(y = y, X = X[, 1:2], Sigma = Sigma, covs = covs, nChunks = 5),
    fastGLS(y = y, X = X, Sigma = Sigma, covs = covs, nChunks = 2)[1:2, ])
})

test_that("genomicControlPValues produces correct new p-values", {
  expect_equal(genomicControlPValues(pVals = .5, nObs = 10)[[1]], 0.5)
  expect_equal(genomicControlPValues(pVals = c(0.25, 0.5), nObs = 10)[[1]],
    c(0.410638105484779, 0.63432328826532))
  expect_equal(genomicControlPValues(pVals = c(0.25, 0.5), nObs = 1e6)[[1]],
    c(0.410588927034021, 0.629472060364479))
})

test_that("genomicControlPValues produces correct inflation factor", {
  expect_equal(genomicControlPValues(pVals = .5, nObs = 10)[[2]], 1)
  expect_equal(genomicControlPValues(pVals = c(0.25, 0.5), nObs = 10)[[2]], 2.04152779518634)
  expect_equal(genomicControlPValues(pVals = c(0.25, 0.5), nObs = 1e6)[[2]], 1.95438310682772)
})

