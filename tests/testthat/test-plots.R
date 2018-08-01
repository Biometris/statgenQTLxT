context("plots")

set.seed(1234)
test_that("qq plot functions properly", {
  pVals <- runif(n = 50)
  p <- qqPlot(pValues = pVals, output = FALSE)
  expect_is(p, "ggplot")
  p1 <-  qqPlot(pValues = pVals, title = "Test title", output = FALSE)
  expect_equal(p1$labels$title, "Test title")
})
