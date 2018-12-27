context("Covariance")

set.seed(1234)

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
               c(3.69773035860965, 4.0458399722506, 3.89885876488689,
                 4.0458399722506, 4.42728292106062, 4.23316194048414,
                 3.89885876488689, 4.23316194048414, 6.02001093783677))
  expect_equal(as.numeric(covunsCov[[2]]),
               c(6.76438475703559, 1.81082720859982, 2.24315513092169,
                 1.81082720859982, 4.74109603467765, 1.41858773775449,
                 2.24315513092169, 1.41858773775449, 3.43705000591464))
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
               c(4.77965289662317, 4.64038160210128, 3.96848926053633,
                 4.64038160210128, 4.59953687220057, 4.3636163363349,
                 3.96848926053633, 4.3636163363349, 6.05946159698509),
               tolerance = 1e-3)
  expect_equal(as.numeric(covpw[[2]]),
               c(5.77371008646564, 1.58047861730721, 2.26327535494345,
                 1.58047861730721, 4.67897152778314, 1.37852313026104,
                 2.26327535494345, 1.37852313026104, 3.40352640707255),
               tolerance = 1e-3)
})

test_that("adding covariates to covPW functions properly", {
  covpwCov <- covPW(Y = Y, K = K, X = X)
  expect_equal(as.numeric(covpwCov[[1]]),
               c(4.06075642932642, 4.23888086262794, 3.85841497799808,
                 4.23888086262794, 4.4701671477277, 4.3593896835162,
                 3.85841497799808, 4.3593896835162, 6.09275913705853),
               tolerance = 1e-3)
  expect_equal(as.numeric(covpwCov[[2]]),
               c(6.29262754840696, 1.70731081261147, 2.28699822432474,
                 1.70731081261147, 4.66326096932245, 1.29752008137592,
                 2.28699822432474, 1.29752008137592, 3.3833683011958),
               tolerance = 1e-3)
})

test_that("option fixDiag functions properly in covPW", {
  expect_warning(covPW(Y = Y, K = K, fixDiag = TRUE),
                 "not implemented yet")
})

test_that("option corMat functions properly in covPW", {
  covPWcor <- covPW(Y = Y, K = K, corMat = TRUE)
  expect_equal(as.numeric(covPWcor[[1]]),
               c(1, -0.718922708200785, -0.210142411208023, -0.718922708200785,
                 1, -0.528493029264564, -0.210142411208023, -0.528493029264564,
                 1), tolerance = 1e-3)
  expect_equal(as.numeric(covPWcor[[2]]),
               c(1, -0.298187737843346, -0.963171857396084, -0.298187737843346,
                 1, 0.543860115826024, -0.963171857396084, 0.543860115826024,
                 1), tolerance = 1e-3)
})

