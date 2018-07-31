context("utils")

test_that("function dfBind copies columns properly", {
  df1 <- data.frame(a = 1:2, b = 1:2)
  df2 <- data.frame(a = 1:2, c = 1:2)
  df3 <- data.frame(c = 1:2, d = 1:2)
  expect_equal(colnames(dfBind(list(df1, df1))), c("a", "b"))
  expect_equal(colnames(dfBind(list(df1, df2))), c("a", "b", "c"))
  expect_equal(colnames(dfBind(list(df1, df3))), c("a", "b", "c", "d"))
  expect_equal(colnames(dfBind(list(df1, df2, df3))), c("a", "b", "c", "d"))
})

test_that("function dfBind inserts NAs for missing columns", {
  df1 <- data.frame(a = 1:2, b = 1:2)
  df2 <- data.frame(a = 1:2, c = 1:2)
  expect_equivalent(unlist(dfBind(list(df1, df2))),
                    c(1, 2, 1, 2, 1, 2, NA, NA, NA, NA, 1, 2))
  expect_equivalent(unlist(dfBind(list(df1, df2, df1))),
                    c(1, 2, 1, 2, 1,2, 1, 2, NA, NA, 1, 2, NA, NA, 1, 2, NA, NA))
})

test_that(paste("function dfBind removes empty data.frames lists from",
                "input before binding"), {
                  df1 <- data.frame(a = 1:2, b = 1:2)
                  expect_equal(dfBind(list(data.frame(), df1)), df1)
                  expect_equal(dfBind(list(df1, data.frame())), df1)
                  expect_null(dfBind(list(data.frame())))
                })

test_that("function matrixRoot functions properly", {
  M1 <- Matrix::Matrix(1:4, nrow = 2)
  M2 <- Matrix::Matrix(c(1:2, 2:1), nrow = 2)
  expect_error(matrixRoot(M1), "should be a symmetric matrix")
  expect_error(matrixRoot(M2), "should be a positive definite matrix")
  expect_equal(as.numeric(matrixRoot(Matrix::crossprod(M2))), c(2, 1, 1, 2))
})


