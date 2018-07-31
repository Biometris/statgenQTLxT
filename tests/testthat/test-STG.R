context("single trait GWAS")

set.seed(1234)
y <- 1:10
X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
Sigma <- matrix(runif(n = 100), nrow = 10)
Sigma <- Sigma %*% t(Sigma)
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

result0 <- runSingleTraitGwas(gData = gDataTest, environments = 1)$GWAResult
result01 <- runSingleTraitGwas(gData = gDataTest)$GWAResult
result1 <- runSingleTraitGwas(gData = gDataTest,
                              environments = 1, covar = "V1")$GWAResult
result2 <- runSingleTraitGwas(gData = gDataTest,
                              environments = 1, snpCov = "M2")$GWAResult
result3 <- runSingleTraitGwas(gData = gDataTest, environments = 1,
                              covar = "V1", snpCov = "M2")$GWAResult

test_that("runSingleTraitGwas produces output with correct dimensions", {
  expect_length(result0, 1)
  expect_length(result01, 2)
})

test_that("runSingleTraitGWas produces correct p-values", {
  expect_equal(result0[[1]]$pValue,
               c(0.517079439679654, 0.91343018536738, 0.628599735847542,
                 0.0807864803940613, 0.857734879152352, 0.0951298087141795,
                 0.609379273189138, 0.999476881270353, 0.41907977041403,
                 0.183886590676029, 0.973491209528092, 0.570647573548852,
                 0.58845656178555, 0.367143146285207, 0.905504194974229))
  expect_equal(result1[[1]]$pValue,
               c(0.269933551571042, 0.984965446392588, 0.648441699544186,
                 0.0626857428326675, 0.866711176390762, 0.120055932703493,
                 0.770457097615474, 0.94062721218024, 0.451231728418542,
                 0.228231902143968, 0.946236384390807, 0.592108077080844,
                 0.61930265357976, 0.405860998685856, 0.911116358438532))
  expect_equal(result2[[1]]$pValue,
               c(0.525193093036643, 0.894587707469969, 0.602679444465947,
                 0.0970424392791122, 0.857734879152351, 0.0928126890021485,
                 0.623720400459772, 0.99078484120111, 0.437029434265121,
                 0.212511818760433, 0.973491209528104, 0.579458353874297,
                 0.648912926295298, 0.36714314628521, 0.918479093445775))
  expect_equal(result3[[1]]$pValue,
               c(0.308728292832839, 0.984965446392575, 0.661002824234323,
                 0.083081633131156, 0.866711176390762, 0.121327025768191,
                 0.804170050729862, 0.927403818290261, 0.465185132885534,
                 0.264951537810541, 0.946236384390817, 0.598229839636104,
                 0.657680037715819, 0.405860998685856, 0.926119151352397))
})
