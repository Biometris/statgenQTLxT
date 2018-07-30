context("Covariance")

## Simulate data.
set.seed(1234)
Vg <- matrix(c(6.06, 4.97, 5.96, 4.97, 4.43, 5.29, 5.96, 5.29, 6.67), nrow = 3)
Ve <- matrix(c(3.28, 0.18, 0.14, 0.18, 4.43, 0.13, 0.14, 0.13, 2.99), nrow = 3)
K <- matrix(c(0.82, 0.007, -0.035, -0.027, -0.035, -0.025, -0.051,
              -0.048, 0.025, -0.032, 0.007, 0.998, 0.025, 0.03, 0.04, 0.022,
              0.212, -0.017, -0.055, -0.011, -0.035, 0.025, 1.124, 0.342, 0.123,
              0.062, 0.097, 0.021, 0.027, -0.047, -0.027, 0.03, 0.342, 1.052,
              0.177, 0.181, 0.004, -0.017, -0.042, -0.036, -0.035, 0.04, 0.123,
              0.177, 1.057, 0.675, 0.016, -0.03, 0.019, -0.014, -0.025, 0.022,
              0.062, 0.181, 0.675, 1.062, -0.019, -0.006, 0.034, -0.009, -0.051,
              0.212, 0.097, 0.004, 0.016, -0.019, 1.136, -0.014, -0.012, -0.037,
              -0.048, -0.017, 0.021, -0.017, -0.03, -0.006, -0.014, 1.115,
              -0.021, -0.076, 0.025, -0.055, 0.027, -0.042, 0.019, 0.034, -0.012,
              -0.021, 1.197, -0.06, -0.032, -0.011, -0.047, -0.036, -0.014,
              -0.009, -0.037, -0.076, -0.06, 0.872), nrow = 10)
S <- matrixRoot(K)
V <- matrixRoot(Vg)
Z <- matrix(rnorm(n = 30), ncol = 3)
Z2 <- matrix(rnorm(n = 30), ncol = 3)
X <- matrix(sample(x = 1:2, size = 10, replace = TRUE), ncol = 1)
Y  <- S %*% Z %*% V + S %*% Z2
rownames(Y) <- rownames(K) <- colnames(K) <- rownames(X) <- paste0("G", 1:10)
colnames(Y) <- paste0("t", 1:3)

covunst <- covUnstr(Y = Y, K = K)
test_that("covUnstr produces correct output structure", {
  expect_is(covunst, "list")
  expect_length(covunst, 2)
  expect_named(covunst, c("Vg", "Ve"))
  expect_is(covunst[[1]], "dpoMatrix")
  expect_is(covunst[[2]], "dpoMatrix")
})

test_that("covUnstr produces correct results", {
  expect_equal(as.numeric(covunst[[1]]),
               c(0.788950280743808, 1.23470182433618, 0.639820421815256,
                 1.23470182433618, 2.06547144502319, 0.31844102893468,
                 0.639820421815256, 0.31844102893468, 4.02050699683009))
  expect_equal(as.numeric(covunst[[2]]),
               c(1.88935541166536, 0.267410948911334, 0.154653069464685,
                 0.267410948911334, 0.572685310452813, 0.431161585873817,
                 0.154653069464685, 0.431161585873817, 0.32584630145632))
})

test_that("adding covariates to covUnstr functions properly", {
  covunsCov <- covUnstr(Y = Y, K = K, X = X)
  expect_equal(as.numeric(covunsCov[[1]]),
               c(0.815888304180514, 1.45137814762539, -0.431898504628126,
                 1.45137814762539, 2.667342335512, -0.513956579069982,
                 -0.431898504628126, -0.513956579069982, 0.985292660263217))
  expect_equal(as.numeric(covunsCov[[2]]),
               c(2.29770467332722, 0.18763949466206, 0.84156097650263,
                 0.18763949466206, 0.357300986190511, 0.374288226470204,
                 0.84156097650263, 0.374288226470204, 0.581257571076995))
})

test_that("option fixDiag functions properly in covUnstr", {
  expect_warning(covUnstr(Y = Y, K = K, fixDiag = TRUE),
                 "not implemented yet")
})

test_that("option veDiag functions properly in covUnstr", {
  covunstVeD <- covUnstr(Y = Y, K = K, VeDiag = TRUE)
  expect_equal(as.numeric(covunstVeD[[1]]),
               c(0.90595431652043, 1.50857194931866, 0.862353505327844,
                 1.50857194931866, 2.55459086889207, 1.04411653929732,
                 0.862353505327844, 1.04411653929732, 4.42906410450837))
  expect_equal(as.numeric(covunstVeD[[2]]),
               c(1.87577258630024, 0, 0, 0, 0.00221657807964868, 0, 0, 0,
                 0.00308230698921089))
})

## Test for covPW give different results in cran check on gitlab then in
## RStudio. Therefore tolerance is added.
covpw <- covPW(Y = Y, K = K)
test_that("covPW produces correct output structure", {
  expect_is(covpw, "list")
  expect_length(covpw, 2)
  expect_named(covpw, c("Vg", "Ve"))
  expect_is(covpw[[1]], "dpoMatrix")
  expect_is(covpw[[2]], "dpoMatrix")
})

test_that("covPW produces correct results", {
  expect_equal(as.numeric(covpw[[1]]),
               c(0.653764571577872, 1.2104145588217, 0.0176273841980339,
                 1.2104145588217, 2.24850717412422, 0.182115998322076,
                 0.0176273841980339, 0.182115998322076, 2.98725317772489),
               tolerance = 1e-5)
  expect_equal(as.numeric(covpw[[2]]),
               c(2.13126590943264, 0.257790587862521, 0.0656276328231413,
                 0.257790587862521, 0.524127476824658, 0.436109022403954,
                 0.0656276328231413, 0.436109022403954, 0.373928447720833),
               tolerance = 1e-5)
})

test_that("adding covariates to covPW functions properly", {
  covpwCov <- covPW(Y = Y, K = K, X = X)
  expect_equal(as.numeric(covpwCov[[1]]),
               c(0.00212922033234732, 0, 0, 0, 2.17712568801086,
                 -0.244324079784721, 0, -0.244324079784721, 1.0843761087518),
               tolerance = 1e-5)
  expect_equal(as.numeric(covpwCov[[2]]),
               c(2.28252977136244, 0.616822992788182, 0.761736077650586,
                 0.616822992788182, 0.349737535988451, 0.337231268115985,
                 0.761736077650586, 0.337231268115985, 0.348508698409624),
               tolerance = 1e-5)
})

test_that("option fixDiag functions properly in covPW", {
  expect_warning(covPW(Y = Y, K = K, fixDiag = TRUE),
                 "not implemented yet")
})

test_that("option corMat functions properly in covPW", {
  covPWcor <- covPW(Y = Y, K = K, corMat = TRUE)
  expect_equal(as.numeric(covPWcor[[1]]),
               c(1, 0.708184590486365, -0.85516690643881, 0.708184590486365,
                 1, -0.971587193712504, -0.85516690643881, -0.971587193712504,
                 1), tolerance = 1e-5)
  expect_equal(as.numeric(covPWcor[[2]]),
               c(1, -0.83472028965713, -0.999997857170303, -0.83472028965713,
                 1, 0.835851093737417, -0.999997857170303, 0.835851093737417,
                 1), tolerance = 1e-5)
})

