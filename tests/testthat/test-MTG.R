context("multi trait GWAS")

set.seed(1234)
pheno <- data.frame(genotype = paste0("G", 1:100), Y)
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
geno <- matrix(sample(x = c(0, 1), size = 300, replace = TRUE), nrow = 100)
colnames(geno) <- rownames(map) <- paste0("M", 1:3)
rownames(geno) <- paste0("G", 1:100)
gDataTest <- createGData(map = map, geno = geno, kin = K,
                         pheno = list(ph1 = pheno, ph2 = pheno),
                         covar = as.data.frame(X))

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
               rep(c(0.00894974481635454, 0.169187547759988, 0.817110552303964),
                   times = 3))
  expect_equal(mtg0$GWAResult[[1]]$effect,
               c(-1.60802146189912, 0.0537719889692854, 0.0673751323310051,
                 -0.967273996067883, -1.06414680197968, -0.409925358018366,
                 -1.68344415285214, -0.0936377753342815, 0.0594676380491173))
  expect_equal(mtg0$GWAResult[[1]]$effectSe,
               c(0.578726031569347, 0.611889320330465, 0.630538683469109,
                 0.538249736599002, 0.571220711368077, 0.589027581808505,
                 0.517023433448998, 0.559488221493538, 0.58226577428052))
})

test_that("option covar functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, covar = "V1")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.0479325785675067, 0.139994166970466, 0.754778566384533),
                   times = 3))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(-1.28701484320345, 0.50679420357873, -0.120718429952228,
                 -0.640348257876892, -0.66244567296074, -0.566378541758748,
                 -1.3801686185568, 0.374915102426538, -0.339987584841663))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(0.60355437361205, 0.608395408459698, 0.626950970194032,
                 0.558520580630795, 0.563085229702065, 0.583341369724312,
                 0.550183761101524, 0.556899668778469, 0.582167242092285))
})

test_that("option snpCov functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, snpCov = "M2")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.0517862491117532, 0.119038956684978, 0.671761680372513),
                   times = 3))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(-1.25324986622524, 0.369619686682055, -0.226185445223459,
                 -0.683423894361723, -0.774536440972262, -0.668159821646221,
                 -1.39774586634064, 0.301170598632323, -0.402804660488232))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(0.608282027316088, 0.608107641339034, 0.630295121522285,
                 0.559208049204777, 0.55890306778119, 0.581877908450072,
                 0.551512715465803, 0.552330037824436, 0.580306949283149))
})

test_that("option subsetmarkers functions properly", {
  expect_error(runMultiTraitGwas(gData = gDataTest, environments = 1,
                                 subsetMarkers = TRUE),
               "markerSubset cannot be empty")
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1,
                           subsetMarkers = TRUE, markerSubset = 1:2)
  expect_equal(nrow(mtg$GWAResult$ph1), 6)
  expect_equivalent(mtg$GWAResult$ph1,
                    mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M3", ])
})

test_that("option MAF functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, MAF = 0.46)
  expect_equal(nrow(mtg$GWAResult$ph1), 6)
  expect_equivalent(mtg$GWAResult$ph1,
                    mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M1", ])
})

test_that("option covModel functions properly", {
  expect_error(runMultiTraitGwas(gData = gDataTest, environments = 1,
                                 fitVarComp = FALSE),
               "Vg should be a matrix")
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, covModel = "pw")
  expect_equal(mtg$GWAResult[[1]]$pValue,
               rep(c(0.0513794317637917, 0.158191909886505, 0.748257521769409),
                   times = 3))
  expect_equal(mtg$GWAResult[[1]]$effect,
               c(-1.42600775669507, 0.235782988067986, -0.0670187725418316,
                 -0.787608968010119, -0.913033845611151, -0.459726378897083,
                 -1.43422632087247, 0.0413914512887092, 0.118624388805461))
  expect_equal(mtg$GWAResult[[1]]$effectSe,
               c(0.600509483021862, 0.624199297148617, 0.631716939323308,
                 0.556763853864026, 0.581134757748847, 0.589479676074537,
                 0.542223018903781, 0.574990785100263, 0.587108570406214))
  mtg1 <- runMultiTraitGwas(gData = gDataTest, environments = 1,
                            covModel = "fa", tolerance = 1)
  expect_equal(mtg1$GWAResult[[1]]$pValue,
               rep(c(0.153967479605108, 0.107136105175992, 0.833495070778579),
                   times = 3))
  expect_equal(mtg1$GWAResult[[1]]$effect,
               c(-1.15122192992848, 0.296248492849938, -0.131864898248064,
                 -0.529505846801976, -0.8765335698888, -0.486917839054459,
                 -1.23854612569167, 0.142975345648095, -0.131740958219656))
  expect_equal(mtg1$GWAResult[[1]]$effectSe,
               c(0.644288853112291, 0.645288084719989, 0.64482614707686,
                 0.591785390081534, 0.599413902726511, 0.599414541454686,
                 0.599129553928024, 0.606593619774449, 0.606666221685391))
})

test_that("option estCom functions properly", {
  mtg <- runMultiTraitGwas(gData = gDataTest, environments = 1, estCom = TRUE)
  expect_equal(mtg$GWAResult$ph1[, 1:9], mtg0$GWAResult$ph1)
  expect_equal(mtg$GWAResult[[1]]$pValCom,
               rep(c(0.00189237221081064, 0.415436068110387, 0.830869076120399),
                   times = 3))
  expect_equal(mtg$GWAResult[[1]]$effsCom,
               rep(c(-1.41536721665423, -0.409874078896011, -0.112816345861181),
                     times = 3))
  expect_equal(mtg$GWAResult[[1]]$effsComSe,
               rep(c(0.439431686764757, 0.481779942001581, 0.505344231337954),
                   times = 3))
  expect_equal(mtg$GWAResult[[1]]$pValQtlE,
               rep(c(0.374723161309463, 0.112615782861975, 0.641637105524885),
                   times = 3))
})
