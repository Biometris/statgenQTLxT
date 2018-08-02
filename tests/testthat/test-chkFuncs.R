context("check Functions")

set.seed(1234)
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- tcrossprod(Sigma)
covs <- matrix(runif(n = 20, max = 100), nrow = 10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(X) <- rownames(Sigma) <- colnames(Sigma) <- rownames(covs) <-
  paste0("G", 1:10)
colnames(X) <- rownames(map) <- paste0("M", 1:3)
gDataTest <- createGData(map = map, geno = X, kin = Sigma,
                         pheno = list(ph1 = pheno, ph2 = pheno),
                         covar = as.data.frame(covs))

test_that("chkGData functions properly", {
  expect_error(chkGData(), "valid gData object")
  expect_error(chkGData(X), "valid gData object")
  expect_silent(chkGData(gDataTest))
  gDataTest2 <- gDataTest
  gDataTest2$markers <- NULL
  expect_error(chkGData(gDataTest2), "should contain markers")
})

test_that("chkNum functions properly", {
  expect_silent(chkNum(1, min = 0, max = 2))
  expect_error(chkNum("a"), "single numerical value")
  expect_error(chkNum(c(1, 2)), "single numerical value")
  expect_error(chkNum(1, min = NULL, max = 0), "smaller than 0")
  expect_error(chkNum(1, min = 2, max = NULL), "greater than 2")
  expect_error(chkNum(1, min = 2, max = 3), "between 2 and 3")
})
