context("EMFA Helper functions")

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

Vg <- matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2)
Ve <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
Dk <- c(1, 2, 3)
vLst <- makeVLst(Vg, Ve, Dk)
test_that("makeVLst functions properly", {
  expect_is(vLst, "list")
  expect_length(vLst, 3)
  expect_is(vLst[[1]], "matrix")
  expect_equal(as.numeric(vLst[[1]]), c(1.6, 0.4, 0.4, 1.6))
  expect_equal(vLst[[2]], 2 * Vg + Ve)
})

vInvLst <- makeVInvLst(Vg, Ve, Dk)
test_that("makeVInvLst functions properly", {
  expect_is(vInvLst, "list")
  expect_length(vInvLst, 3)
  expect_is(vInvLst[[1]], "matrix")
  expect_equal(as.numeric(vInvLst[[1]]),
               c(0.666666666666667, -0.166666666666667, -0.166666666666667,
                 0.666666666666667))
  expect_equal(vInvLst[[2]], solve(2 * Vg + Ve))
})
