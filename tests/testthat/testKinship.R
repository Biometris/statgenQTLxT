context("Test of kinship functions")

test_that("kinship functions give correct output", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(as.matrix(astle(X)), matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))
  expect_equivalent(as.matrix(GRM(X)), matrix(c(0.5, -0.5, -0.5, 0.5), nrow = 2))
  expect_equivalent(as.matrix(IBS(X)), matrix(c(1, 0, 0, 1), nrow = 2))
  expect_equivalent(as.matrix(vanRaden(X)), matrix(c(2/3, -2/3, -2/3, 2/3), nrow = 2))
})

test_that("kinship functions give correct output with user definded denominator", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_equivalent(as.matrix(astle(X, denominator = 4)), matrix(c(1/3, -1/3, -1/3, 1/3), nrow = 2))
  expect_equivalent(as.matrix(GRM(X, denominator = 4)), matrix(c(0.25, -0.25, -0.25, 0.25), nrow = 2))
  expect_equivalent(as.matrix(IBS(X, denominator = 4)), matrix(c(0.5, 0, 0, 0.5), nrow = 2))
  expect_equivalent(as.matrix(vanRaden(X, denominator = 4)), matrix(c(1/8, - 1/8, -1/8, 1/8), nrow = 2))
})
