context("kinship functions")

test_that("kinship functions give correct output", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X), matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))
  expect_equivalent(GRMCPP(X), matrix(c(0.5, -0.5, -0.5, 0.5), nrow = 2))
  expect_equivalent(IBSCPP(X), matrix(c(1, 0, 0, 1), nrow = 2))
  expect_equivalent(vanRadenCPP(X), matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))
})

test_that(paste("kinship functions give correct output with user",
                "definded denominator"), {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(astleCPP(X, denom = 4),
                    matrix(c(1/3, -1/3, -1/3, 1/3), nrow = 2))
  expect_equivalent(GRMCPP(X, denom = 4),
                    matrix(c(0.25, -0.25, -0.25, 0.25), nrow = 2))
  expect_equivalent(IBSCPP(X, denom = 4),
                    matrix(c(0.5, 0, 0, 0.5), nrow = 2))
  expect_equivalent(vanRadenCPP(X, denom = 4),
                    matrix(c(1/8, -1/8, -1/8, 1/8), nrow = 2))
})

set.seed(1234)
map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                  pos = c(1, 4, 17, 2, 3, 5))
posCor <- c(3, 8, 13, 1, 1.5, 2)
X <- matrix(runif(n = 180), nrow = 60)
X <- X / rowSums(X)

dim(X) <- c(10, 6, 3)
test_that("multiAllKinCPP functions properly", {
  K <- multiAllKinCPP(x = X, posCor = posCor)
  expect_is(K, "matrix")
  expect_equal(dim(K), c(10, 10))
  expect_equal(diag(K), rep(x = 1, times = 10))
  expect_equal(K[, 1],
               c(1, 0.318199140315202, 0.353188959968734, 0.356389520929043,
                 0.365703566465929, 0.344204968497873, 0.363187385587088,
                 0.353630750320603, 0.303310956779327, 0.383196299221719))
  K2 <- multiAllKinCPP(x = X, posCor = posCor, denom = 1)
  expect_equal(K2 / 28.5, K)
})
