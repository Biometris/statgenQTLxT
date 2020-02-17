load(file = "testdata.rda")

### Test multi trait GWAS.

set.seed(1234)
pheno <- data.frame(genotype = rownames(Y), Y)
map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
geno <- matrix(sample(x = c(0, 1), size = 300, replace = TRUE), nrow = 100)
colnames(geno) <- rownames(map) <- paste0("M", 1:3)
rownames(geno) <- rownames(Y)
gDataTest <- createGData(map = map, geno = geno, kin = K,
                         pheno = list(ph1 = pheno, ph2 = pheno),
                         covar = as.data.frame(W))

# Check output structure.

expect_error(runMultiTraitGwas(gData = gDataTest), "Trial cannot be NULL")

mtg0 <- runMultiTraitGwas(gData = gDataTest, trials = 1)

expect_true(inherits(mtg0, "GWAS"))
expect_equal(length(mtg0), 5)
expect_equal(names(mtg0),
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_true(inherits(mtg0$GWAResult, "list"))
expect_equal(length(mtg0$GWAResult), 1)
expect_equal(names(mtg0$GWAResult), "ph1")

# Check results.

expect_equal(mtg0$GWAResult[[1]]$pValue,
             rep(c(0.564371004965954, 0.940801610756639, 0.332245567117611),
                 times = 3))
expect_equal(mtg0$GWAResult[[1]]$effect,
             c(0.0167297922953778, -0.0149118755656098, 0.121523992392644,
               -0.24418433632583, -0.0684450538191805, -0.0738500285513413,
               -0.0445236264165599, 0.0466707878304148, 0.182116758308792))
expect_equal(mtg0$GWAResult[[1]]$effectSe,
             c(0.141486409847841, 0.140273821295455, 0.14149047684672,
               0.197714552966762, 0.196078915075951, 0.198247455496308,
               0.15023431551672, 0.148812963838418, 0.149778360788749))

# Add covariates.

mtg1 <- runMultiTraitGwas(gData = gDataTest, trials = 1, covar = "V1")
expect_equal(mtg1$GWAResult[[1]]$pValue,
             rep(c(1, 0.960604256923414, 0.318227779346518), times = 3))
expect_equal(mtg1$GWAResult[[1]]$effect,
             c(-0.0400403145641162, -0.0140461915962329, 0.123883152061841,
               -2.50466324253816, -0.0347502741529288, -0.0751942786184463,
               -1.26893671653622, 0.0565051030727486, 0.184737941923979))
expect_equal(mtg1$GWAResult[[1]]$effectSe,
             c(NaN, 0.142612687595557, 0.14230199541243, 1438659.83234064,
               0.197967003080701, 0.197978147592246, 1233485.45212262,
               0.150877382209886, 0.150493360202025))

# Add SNP covariates.

mtg2 <- runMultiTraitGwas(gData = gDataTest, trials = 1, snpCov = "M2")
expect_equal(mtg2$GWAResult[[1]]$pValue,
             rep(c(0.591629057264197, 0.929162416087565, 0.381571064212693),
                 times = 3))
expect_equal(mtg2$GWAResult[[1]]$effect,
             c(0.0188979667278949, -0.0149368688983746, 0.122019078604419,
               -0.240151777405656, -0.0753705460378873, -0.0620656988090511,
               -0.0512324342564622, 0.0480887401914931, 0.176639092147782))
expect_equal(mtg2$GWAResult[[1]]$effectSe,
             c(0.143697512244861, 0.14093370713494, 0.142512013201738,
               0.201709601438192, 0.197821205742041, 0.200265131684195,
               0.152406568586504, 0.149552893298575, 0.151065022278834))

# Test option subsetMarkers

expect_error(runMultiTraitGwas(gData = gDataTest, trials = 1,
                               subsetMarkers = TRUE),
             "markerSubset cannot be empty")
mtg3 <- runMultiTraitGwas(gData = gDataTest, trials = 1,
                          subsetMarkers = TRUE, markerSubset = 1:2)
expect_equal(nrow(mtg3$GWAResult$ph1), 6)
expect_equivalent(mtg3$GWAResult$ph1,
                  mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M3", ])

# Test option MAF.

mtg4 <- runMultiTraitGwas(gData = gDataTest, trials = 1, MAF = 0.45)
expect_equal(nrow(mtg4$GWAResult$ph1), 6)
expect_equivalent(mtg4$GWAResult$ph1,
                  mtg0$GWAResult$ph1[mtg0$GWAResult$ph1$snp != "M1", ])

# Test option covModel.

expect_error(runMultiTraitGwas(gData = gDataTest, trials = 1,
                               fitVarComp = FALSE),
             "Vg should be a matrix")
mtg5 <- runMultiTraitGwas(gData = gDataTest, trials = 1, covModel = "pw")
expect_equal(mtg5$GWAResult[[1]]$pValue,
             rep(c(0.584193379521886, 0.917792804538605, 0.215174046880372),
                 times = 3))
expect_equal(mtg5$GWAResult[[1]]$effect,
             c(0.0321190670026913, -0.0276679839039544, 0.127800258683102,
               -0.223755785437776, -0.105539529188826, -0.0820764038081917,
               -0.0366464885798382, 0.0205386668807394, 0.21330502269409))
expect_equal(mtg5$GWAResult[[1]]$effectSe,
             c(0.143190492553253, 0.141877560601682, 0.142670002350013,
               0.197681879042029, 0.196006499919309, 0.19802527279839,
               0.153014644408791, 0.151512670668318, 0.151760600994316))

mtg6 <- runMultiTraitGwas(gData = gDataTest, trials = 1, covModel = "fa",
                          tolerance = 1)
expect_equal(mtg6$GWAResult[[1]]$pValue,
             rep(c(0.62520660892386, 0.893120546011254, 0.220861887600226),
                 times = 3))
expect_equal(mtg6$GWAResult[[1]]$effect,
             c(0.0264111446640612, -0.0422680217968601, 0.138984437262161,
               -0.217627065471035, -0.107145409801779, -0.061239108238967,
               -0.0340881382484514, 0.0332767984091363, 0.217851758849069))
expect_equal(mtg6$GWAResult[[1]]$effectSe,
             c(0.144131560902014, 0.142756752545639, 0.14323134670053,
               0.202833551209507, 0.200907212892662, 0.201658772009655,
               0.151874982098179, 0.150403664754337, 0.150707331101958))

# Test option estCom.

mtg7 <- runMultiTraitGwas(gData = gDataTest, trials = 1, estCom = TRUE)
expect_equal(mtg7$GWAResult$ph1[, 1:9], mtg0$GWAResult$ph1)
expect_equal(mtg7$GWAResult[[1]]$pValCom,
             rep(c(0.850460837572108, 0.935381381322014, 0.23123256007111),
                 times = 3))
expect_equal(mtg7$GWAResult[[1]]$effsCom,
             rep(c(-0.0220950879871188, 0.00942387367522931, 0.140597250775363),
                 times = 3))
expect_equal(mtg7$GWAResult[[1]]$effsComSe,
             rep(c(0.116903742175511, 0.115941743450033, 0.117284736362518),
                 times = 3))
expect_equal(mtg7$GWAResult[[1]]$pValQtlE,
             rep(c(0.367863250786882, 0.822802895922502, 0.371379782602173),
                 times = 3))
