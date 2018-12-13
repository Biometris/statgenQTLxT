context("Covariance")

covunst <- covUnstr(Y = Y, K = K)
test_that("covUnstr produces correct output structure", {
  expect_is(covunst, "list")
  expect_length(covunst, 2)
  expect_named(covunst, c("Vg", "Ve"))
  expect_is(covunst[[1]], "matrix")
  expect_is(covunst[[2]], "matrix")
})

test_that("covUnstr produces correct results", {
  expect_equal(as.numeric(covunst[[1]]),
               c(4.2913219214206, 4.41303217609796, 3.98917611215851,
                 4.41303217609796, 4.54150675723924, 4.18905792653371,
                 3.98917611215851, 4.18905792653371, 5.97996073124397))
  expect_equal(as.numeric(covunst[[2]]),
               c(6.3354652747199, 1.68406706569666, 2.24476139489599,
                 1.68406706569666, 4.80838701002544, 1.52825888742586,
                 2.24476139489599, 1.52825888742586, 3.46205962753695))
})

test_that("adding covariates to covUnstr functions properly", {
  covunsCov <- covUnstr(Y = Y, K = K, X = X)
  expect_equal(as.numeric(covunsCov[[1]]),
               c(4.2913219214203, 4.41303217609752, 3.98917611215893,
                 4.41303217609752, 4.54150675723864, 4.1890579265338,
                 3.98917611215893, 4.1890579265338, 5.97996073124403))
  expect_equal(as.numeric(covunsCov[[2]]),
               c(6.33546527471985, 1.6840670656963, 2.24476139489616,
                 1.6840670656963, 4.80838701002509, 1.52825888742582,
                 2.24476139489616, 1.52825888742582, 3.46205962753691))
})

test_that("option fixDiag functions properly in covUnstr", {
  expect_warning(covUnstr(Y = Y, K = K, fixDiag = TRUE),
                 "not implemented yet")
})

test_that("option veDiag functions properly in covUnstr", {
  covunstVeD <- covUnstr(Y = Y, K = K, VeDiag = TRUE)
  expect_equal(as.numeric(covunstVeD[[1]]),
               c(6.79924232675014, 6.48978847589632, 6.651961918426,
                 6.48978847589632, 6.21530475446956, 6.19220313491186,
                 6.651961918426, 6.19220313491186, 7.68818904405218))
  expect_equal(as.numeric(covunstVeD[[2]]),
               c(4.44472679015, 0, 0, 0, 3.62606231880493, 0, 0, 0,
                 2.25759456110549))
})

## Test for covPW give different results in cran check on gitlab compared to
## RStudio. Therefore tolerance is added.
covpw <- covPW(Y = Y, K = K)
test_that("covPW produces correct output structure", {
  expect_is(covpw, "list")
  expect_length(covpw, 2)
  expect_named(covpw, c("Vg", "Ve"))
  expect_is(covpw[[1]], "matrix")
  expect_is(covpw[[2]], "matrix")
})

test_that("covPW produces correct results", {
  expect_equal(as.numeric(covpw[[1]]),
               c(4.77874019368121, 4.64063366081611, 3.98034434230101,
                 4.64063366081611, 4.5977520787454, 4.36565379102361,
                 3.98034434230101, 4.36565379102361, 6.05931461904336),
               tolerance = 1e-3)
  expect_equal(as.numeric(covpw[[2]]),
               c(5.77371008646564, 1.58138500737496, 2.25287028616509,
                 1.58138500737496, 4.67897152778314, 1.37706849127846,
                 2.25287028616509, 1.37706849127846, 3.40352640707255),
               tolerance = 1e-3)
})

test_that("adding covariates to covPW functions properly", {
  covpwCov <- covPW(Y = Y, K = K, X = X)
  expect_equal(as.numeric(covpwCov[[1]]),
               c(4.77874019368084, 4.64063366081565, 3.98034434230131,
                 4.64063366081565, 4.59775207874487, 4.36565379102375,
                 3.98034434230131, 4.36565379102375, 6.05931461904335),
               tolerance = 1e-3)
  expect_equal(as.numeric(covpwCov[[2]]),
               c(5.77371008646564, 1.58138500737467, 2.25287028616532,
                 1.58138500737467, 4.67897152778313, 1.37706849127842,
                 2.25287028616532, 1.37706849127842, 3.40352640707255),
               tolerance = 1e-3)
})

test_that("option fixDiag functions properly in covPW", {
  expect_warning(covPW(Y = Y, K = K, fixDiag = TRUE),
                 "not implemented yet")
})

test_that("option corMat functions properly in covPW", {
  covPWcor <- covPW(Y = Y, K = K, corMat = TRUE)
  expect_equal(as.numeric(covPWcor[[1]]),
               c(1, -0.716293991666446, -0.215541722334213, -0.716293991666446,
                 1, -0.527005339751727, -0.215541722334213, -0.527005339751727,
                 1), tolerance = 1e-3)
  expect_equal(as.numeric(covPWcor[[2]]),
               c(1, -0.307156052294611, -0.964460926783025, -0.307156052294611,
                 1, 0.547692563945474, -0.964460926783025, 0.547692563945474,
                 1), tolerance = 1e-3)
})

