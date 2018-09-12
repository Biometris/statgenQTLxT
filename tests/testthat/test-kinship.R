context("kinship functions")

test_that("kinship functions give correct output", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X),c(2/3, -2/3, -2/3, 2/3))
  expect_equivalent(GRMCPP(X), c(0.5, -0.5, -0.5, 0.5))
  expect_equivalent(IBSCPP(X), c(1, 0, 0, 1))
  expect_equivalent(vanRadenCPP(X), c(2/3, -2/3, -2/3, 2/3))
})

test_that(paste("kinship functions give correct output with user",
                "definded denominator"), {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X, denom = 4), c(1/3, -1/3, -1/3, 1/3))
  expect_equivalent(GRMCPP(X, denom = 4), c(0.25, -0.25, -0.25, 0.25))
  expect_equivalent(IBSCPP(X, denom = 4), c(0.5, 0, 0, 0.5))
  expect_equivalent(vanRadenCPP(X, denom = 4), c(1/8, -1/8, -1/8, 1/8))
})

test_that("function kinship functions properly for 2d markers", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equal(kinship(X = X, method = "astle"), astleCPP(X))
  expect_equal(kinship(X = X, method = "GRM"), GRMCPP(X))
  expect_equal(kinship(X = X, method = "IBS", denominator = 2),
               IBSCPP(X, denom = 2))
  expect_equal(kinship(X = X, method = "vanRaden", denominator = 3),
               vanRadenCPP(X, denom = 3))
})

set.seed(1234)
map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                  pos = c(1, 4, 17, 2, 3, 5))
posCor <- c(3, 8, 13, 1, 1.5, 2)
X <- matrix(runif(n = 180), nrow = 60)
X <- X / rowSums(X)
dim(X) <- c(10, 6, 3)
K0 <- multiAllKinCPP(x = X, posCor = posCor)
test_that("multiAllKinCPP functions properly", {
  expect_is(K0, "matrix")
  expect_equal(dim(K0), c(10, 10))
  expect_equal(diag(K0), rep(x = 1, times = 10))
  expect_equal(K0[, 1],
               c(1, 0.318199140315202, 0.353188959968734, 0.356389520929043,
                 0.365703566465929, 0.344204968497873, 0.363187385587088,
                 0.353630750320603, 0.303310956779327, 0.383196299221719))
  K1 <- multiAllKinCPP(x = X, posCor = posCor, denom = 1)
  expect_equal(K1 / 28.5, K0)
})

test_that("function kinship functions properly for 3d markers", {
  K2 <- kinship(X = X, map = map, method = "multiAllKin")
  expect_equal(K2, K0)
  K3 <- kinship(X = X, map = map, method = "multiAllKin", denominator = 1)
  expect_equal(K3 / 28.5, K2)
})

