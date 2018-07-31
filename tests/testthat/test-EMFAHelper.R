context("EMFA Helper functions")

test_that("startValPW functions properly", {
  C <- Matrix::Matrix(c(2,1,1,2), nrow = 2)
  expect_null(startValPW(C = C, m = 0, p = 2)[[1]])
  startVal <- startValPW(C = C, m = 1, p = 2)
  expect_is(startVal, "list")
  expect_named(startVal, c("W", "P"))
  expect_is(startVal[[1]], "Matrix")
  expect_equal(dim(startVal[[1]]), c(2, 1))
  expect_equal(dim(startVal[[2]]), c(2, 2))
  expect_equal(as.numeric(startVal[[1]]),
               c(-0.577350269189626, 0.577350269189626))
  expect_equal(as.numeric(startVal[[2]]), c(3, 0, 0, 3))
})

C <- Matrix::Matrix(c(2,1,1,2), nrow = 2)
Omega <- Matrix::Matrix(c(0.5, 0.3, 0.3, 0.5), nrow = 2)
startVal <- startValPW(C = C, m = 1, p = 2)
W <- startVal$W
P <- startVal$P
pr0 <- updatePrec(m = 0, p = 2, Omega = Omega, W = W, P = P, het = FALSE,
                  maxDiag = 1e5)
pr1 <- updatePrec(m = 1, p = 2, Omega = Omega, W = W, P = P, het = FALSE,
                  maxDiag = 1e5)
test_that("updatePrec produces correct output structure", {
  expect_is(pr1, "list")
  expect_named(pr1, c("CNew", "WNew", "PNew"))
  expect_is(pr1[[1]], "dsyMatrix")
  expect_is(pr1[[2]], "dgeMatrix")
  expect_is(pr1[[3]], "ddiMatrix")
  expect_equal(dim(pr1[[1]]), c(2, 2))
  expect_equal(dim(pr1[[2]]), c(2, 1))
  expect_equal(dim(pr1[[3]]), c(2, 2))
  expect_null(pr0[[2]])
  expect_equal(pr0[[1]], pr0[[3]])
})

test_that("updatePrec functions properly", {
  expect_equal(as.numeric(pr0[[1]]), c(2, 0, 0, 2))
  expect_equal(as.numeric(pr1[[1]]), c(3.125, -1.875, -1.875, 3.125))
  expect_equal(as.numeric(pr1[[2]]), c(0.547722557505166, 0.547722557505166))
  expect_equal(as.numeric(pr1[[3]]), c(5, 0, 0, 5))
})

test_that("option het in updatePrec functions properly", {
  pr0a <- updatePrec(m = 0, p = 2, Omega = Omega, W = W, P = P, het = TRUE,
                     maxDiag = 1e5)
  pr1a <- updatePrec(m = 1, p = 2, Omega = Omega, W = W, P = P, het = TRUE,
                     maxDiag = 1e5)
  expect_equal(as.numeric(pr0a[[1]]), c(2, 0, 0, 2))
  expect_equal(as.numeric(pr1a[[1]]),
               c(2, 2.57226643639495e-09, 2.57226643639495e-09, 2))
  expect_equal(as.numeric(pr1a[[2]]),
               c(-2.53587580354153e-05, 2.53587580354153e-05))
  expect_equal(as.numeric(pr1a[[3]]),
               c(2.00000000257227, 0, 0, 2.00000000257227))
})

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

Vg2 <- Matrix::Matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2)
Ve2 <- Matrix::Matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
vLst2 <- makeVLst(Vg2, Ve2, Dk)
vInvLst2 <- makeVInvLst(Vg2, Ve2, Dk)
test_that("makeVLst and makeVInvLst retain class", {
  expect_is(vLst2[[1]], "Matrix")
  expect_is(vInvLst2[[1]], "Matrix")
})

Y <- Matrix::Matrix(1:6, nrow = 2)
X <- Matrix::Matrix(1:3, nrow = 1)
test_that("LLDiag functions properly", {
  llDiag <- LLDiag(Y = Y, vLst = vLst2, vInvLst = vInvLst2)
  expect_is(llDiag, "numeric")
  expect_length(llDiag, 1)
  expect_equal(llDiag, -15.5221797651405)
  llDiag2 <- LLDiag(Y = Y, X = X, vLst = vLst2, vInvLst = vInvLst2)
  expect_equal(llDiag2, -2.38746869986328)
})

test_that("LLQuadFormDiag functions properly", {
  llQFDiag <- LLQuadFormDiag(Y = Y, vInvLst = vInvLst2)
  expect_is(llQFDiag, "numeric")
  expect_length(llQFDiag, 1)
  expect_equal(llQFDiag, 26.5208333333333)
  llQFDiag2 <- LLQuadFormDiag(Y = Y, X = X, vInvLst = vInvLst2)
  expect_equal(llQFDiag2, -663.561705739926)
})

