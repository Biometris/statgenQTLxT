# context("EMFA Helper functions")
#
# Dk <- c(1, 2, 3)
# Vg2 <- Matrix::Matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2)
# Ve2 <- Matrix::Matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
# vInvLst2 <- makeVInvLst(Vg2, Ve2, Dk)
#
# vInvArr <- array(sapply(vInvLst2, as.matrix), dim = c(2, 2, 3))
# test_that("LLQuadFormDiag functions properly", {
#   llQFDiag <- LLQuadFormDiagCPP(y = as.matrix(Y), vInv = vInvArr)
#   expect_is(llQFDiag, "numeric")
#   expect_length(llQFDiag, 1)
#   expect_equal(llQFDiag, 26.5208333333333)
#   llQFDiag2 <- LLQuadFormDiagCPP(y = as.matrix(Y), vInv = vInvArr,
#                                  size_param = as.matrix(X))
#   expect_equal(llQFDiag2, -663.561705739926)
# })
#
# vInvArr2 <- vInvArr[1, 1, , drop = FALSE]
# Y2 <- as.matrix(Y[1, , drop = FALSE])
# test_that("LLQuadFormDiag functions properly for dimension 1 by 1", {
#   llQFDiag <- LLQuadFormDiagCPP(y = Y2, vInv = vInvArr2)
#   expect_equal(llQFDiag, 14.3541666666667)
#   llQFDiag2 <- LLQuadFormDiagCPP(y = Y2, vInv = vInvArr2,
#                                  size_param = as.matrix(X))
#   expect_equal(llQFDiag2, -184.528147600263)
# })
#
