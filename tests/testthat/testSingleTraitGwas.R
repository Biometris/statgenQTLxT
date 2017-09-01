context("Test single trait Gwas and helper functions.")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)

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

