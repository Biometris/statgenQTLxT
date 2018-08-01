context("multi trait GWAS")

set.seed(1234)
y <- 1:10
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

mtg0 <- runMultiTraitGwas(gData = gDataTest, environments = 1)
test_that("runMultiTraitGwas produces correct output structure", {
  expect_error(runMultiTraitGwas(gData = gDataTest),
               "Environment cannot be NULL")
  expect_is(mtg0, "GWAS")
  expect_length(mtg0, 5)
  expect_named(mtg0, c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
  expect_is(mtg0$GWAResult, "list")
  expect_length(mtg0$GWAResult, 1)
  expect_named(mtg0$GWAResult, "ph1")
})

test_that("runMultiTraitGwas produces correct values", {
  expect_equal(mtg0$GWAResult[[1]]$pValue,
               rep(c(0.264438604699095, 0.760777732822648, 0.713811311590879),
                   times = 5))
  expect_equal(mtg0$GWAResult[[1]]$effect,
               c(0.651462325998874, -0.0963157332445604, -0.45493578848263,
                 -2.85786967095735, -0.162111751230252, 1.61881209347749,
                 -0.477917062173166, -0.570596690857517, -1.23364774969637,
                 0.553318859071792, -0.255894304568591, -0.505407113745332,
                 -2.44953282738719, -2.01618256813099, -0.252750203906081))
  expect_equal(mtg0$GWAResult[[1]]$effectSe,
               c(1.19818120578017, 1.12056339432964, 1.12059264905764,
                 2.50600541312395, 2.23502213484021, 2.29185399382522,
                 2.41738666691005, 2.07090388480029, 2.17159964975999,
                 1.39841894683381, 1.08911387624451, 1.19340270892044,
                 2.49754799236346, 1.92456105970911, 2.12845105537178))
})

test_that("option covar functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, covar = "V1")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.229590455157771, 0.970705786567212, 0.632754027316857),
                   times = 5))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(1.24407082741394, -0.0191668919515887, -0.449811400736985,
                 -2.11690557400844, 0.116036019105379, 1.47774867791872,
                 0.872789917032494, -0.278696370407305, -1.47736503418324,
                 1.52752982794928, -0.157599471603168, -0.794068954030965,
                 -0.410650887299702, -1.5005581682661, -0.32042750116859))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(1.37284773972262, 1.1614288682791, 1.155193590822,
                 2.89561066366456, 2.59137053168845, 2.44475192347802,
                 2.60178657731, 2.46554187760066, 2.19893907635258,
                 1.40795388511472, 1.41713017419629, 1.17847756124152,
                 2.47435112046887, 2.56797495765052, 2.08584944982048))
})

test_that("option snpCov functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, snpCov = "M2")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.390330854695164, 0.958052185333562, 0.618568268857209),
                   times = 5))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(0.666921304210177, -0.096761685073045, -0.511353226667357,
                 -1.82936460662647, 0.167145356727907, 1.80874609536277,
                 0.886097306865886, -0.152421418895658, -1.37179452867082,
                 1.34760341053236, -0.0565518032016152, -0.7238188354907,
                 -0.443431954426977, -1.42206386837535, -0.380927292648401))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(1.2422052751751, 1.15740884785855, 1.19516967550706,
                 2.61780938977131, 2.52858093680548, 2.50614744034038,
                 2.35017324374819, 2.35410574979857, 2.23339299935478,
                 1.27148973420054, 1.3382296126582, 1.17111527315,
                 2.19354619725693, 2.36864757472567, 2.02663639387221))
})

test_that("option subsetmarkers functions properly", {
  expect_error(runMultiTraitGwas(gData = gDataTest, environments = 1,
                                 subsetMarkers = TRUE),
               "markerSubset cannot be empty")
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1,
                           subsetMarkers = TRUE, markerSubset = 1:2)
  expect_equal(nrow(mtg$GWAResult$ph1), 10)
  expect_equivalent(mtg$GWAResult$ph1,
                    mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M3", ])
})

test_that("option MAF functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, MAF = 0.4)
  expect_equal(nrow(mtg$GWAResult$ph1), 10)
  expect_equivalent(mtg$GWAResult$ph1,
                    mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M1", ])
})

test_that("option covModel functions properly", {
  expect_error(runMultiTraitGwas(gData = gDataTest, environments = 1,
                                 fitVarComp = FALSE),
               "Vg should be a matrix")
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, covModel = "pw")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.571618685792894, 0.934957202566941, 0.799194594386479),
                   times = 5))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(0.560441865812252, -0.202006679671093, -0.460510750784809,
                 -2.17190695243556, -0.0364055857309747, 1.78835682077619,
                 0.433723230352175, -0.344132208412084, -1.0481415913455,
                 1.1338358674987, -0.0286402010583331, -0.474527359983557,
                 -0.850317469397694, -1.29719883559961, 0.175040316689073))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(1.04280266783828, 0.911605531523456, 0.947834846772988,
                 1.47872735314371, 1.1831429265397, 1.30607970264501,
                 1.55892899510557, 1.3046535188873, 1.39662404885936,
                 0.82620078752817, 0.772582281121678, 0.77261489702717,
                 1.44446191547629, 1.35071941393981, 1.35077627811758))
  mtg1 <- runMultiTraitGwas(gData = gDataTest, environments = 1,
                            covModel = "fa", tolerance = 1)
  expect_equal(mtg1$GWAResult[[1]]$pValue,
               rep(c(0.737435444098133, 0.924722166075568, 0.867268991529069),
                   times = 5))
  expect_equal(mtg1$GWAResult[[1]]$effect,
               c(-0.37107888621318, -0.112591521081918, 0.234731362844117,
                 -0.920420904286224, -0.0241792058948211, 1.01546383054652,
                 -0.360996718362455, -0.126860389072551, -0.407499160652376,
                 0.0799843969241693, -0.0354023233107799, 0.239237527985091,
                 -1.68524110911562, -1.30044842458866, 0.725760186454197))
  expect_equal(mtg1$GWAResult[[1]]$effectSe,
               c(0.7697760390041, 0.586409983013696, 0.669688205356899,
                 0.928809299358635, 0.764521762981082, 0.824386567914816,
                 1.27382738346414, 1.12652923886213, 1.16426330694102,
                 0.785758003251032, 0.598938105720756, 0.684310655259404,
                 1.45316810944861, 1.31643360203728, 1.3412471291551))
})

test_that("option estCom functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, estCom = TRUE)
  expect_equal(mtg$GWAResult$ph1[1:9], mtg0$GWAResult$ph1)
  expect_equal(mtg$GWAResult[[1]]$pValCom,
               rep(c(0.753742167795203, 0.411338088691726, 0.464498438438716),
                   times = 5))
  expect_equal(mtg$GWAResult[[1]]$effsCom,
               rep(c(0.1779474955898, -0.376742575564495, -0.366029632818305),
                     times = 5))
  expect_equal(mtg$GWAResult[[1]]$effsComSe,
               rep(c(0.863716901921067, 0.700585528254365, 0.763875628958801),
                   times = 5))
  expect_equal(mtg$GWAResult[[1]]$pValQtlE,
               rep(c(0.178732679052409, 0.747264697590076, 0.667677655472902),
                   times = 5))
})
