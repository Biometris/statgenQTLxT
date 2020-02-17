set.seed(1234)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- tcrossprod(Sigma)
rownames(Sigma) <- colnames(Sigma) <- paste0("G", 1:10)
pheno <- data.frame(genotype = paste0("G", 1:10),
                    matrix(rnorm(50, mean = 10, sd = 2), nrow = 10))
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
rownames(map) <- paste0("M", 1:3)
XArr <- matrix(runif(n = 90), nrow = 30)
XArr <- XArr / rowSums(XArr)
dim(XArr) <- c(10, 3, 3)
rownames(XArr) <- paste0("G", 1:10)
colnames(XArr) <- paste0("M", 1:3)
dimnames(XArr)[[3]] <- paste0("A", 1:3)
gDataTestIBD <- createGData(map = map, geno = XArr, kin = Sigma, pheno = pheno)

## Test matrix plot

stgIBD <- runSingleTraitGwasIBD(gDataTestIBD, traits = "X1")
p <- plot(stgIBD, plotType = "matrix", output = FALSE)
expect_true(inherits(p, "ggplot"))
p1 <- plot(stgIBD, plotType = "matrix", xLab = "labx", yLab = "laby",
           output = FALSE)
expect_equal(p1$labels$x, "labx")
expect_equal(p1$labels$y, "laby")

