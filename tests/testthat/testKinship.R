X <- matrix(c(1, 0, 0, 1), nrow = 2)

test_that("Test output kinship functions.", {
  expect_equal(astle(X), matrix(c(2 / 3, -2 / 3, -2 / 3, 2 / 3), nrow = 2))
  expect_equal(GRM(X), matrix(c(0.5, 0, 0, 0.5), nrow = 2))
  expect_equal(IBS(X), matrix(c(1, 0, 0, 1), nrow = 2))
  expect_equal(vanRaden(X), matrix(c(2 / 3, -2 / 3, -2 / 3, 2 / 3), nrow = 2))
})

