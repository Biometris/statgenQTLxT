context("FA Helper functions")

test_that("vecInvDiag functions properly", {
  vid <- vecInvDiag(1:2, 3:5)
  expect_is(vid, "matrix")
  expect_equal(dim(vid), c(3, 2))
  expect_equal(as.numeric(vid),
               c(0.25, 0.2, 0.166666666666667, 0.142857142857143,
                 0.111111111111111, 0.0909090909090909))
})

test_that("tracePInvDiag functions properly", {
  tpid <- tracePInvDiag(1:2, 3:5)
  expect_is(tpid, "numeric")
  expect_length(tpid, 2)
  expect_equal(tpid, c(0.616666666666667, 0.344877344877345))
})

