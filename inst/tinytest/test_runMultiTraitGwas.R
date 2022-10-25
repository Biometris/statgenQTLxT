load(file = "testdata.rda")

gDataTest <- statgenGWAS::createGData(geno = X, map = map, pheno = Y, covar = W)

### General input checks
expect_error(runMultiTraitGwas(1),
             "gData should be a valid gData object")
expect_error(runMultiTraitGwas(gDataTest, trials = "a"),
             "trials should be in pheno")


### Test runMultiTraitGwas with single kinship.

mtg0 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"))

result1 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             covModel = "pw", covar = "V1")[["GWAResult"]]
result1a <- runMultiTraitGwas(gData = gDataTest, traits = 2:3,
                              covModel = "pw", covar = 1)[["GWAResult"]]
result1b <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              covModel = "pw", covar = 1)[["GWAResult"]]
result2 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             covModel = "fa", snpCov = "M002")[["GWAResult"]]
result3 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             covModel = "fa", covar = "V1",
                             snpCov = "M002")[["GWAResult"]]

## Check output structure.

expect_inherits(mtg0, "GWAS")
expect_equal(length(mtg0), 5)
expect_equal(names(mtg0),
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal(names(mtg0[["GWASInfo"]]),
             c("call", "MAF", "thrType", "GLSMethod", "covModel",
               "varComp", "genomicControl", "inflationFactor"))
expect_equal(length(mtg0[["GWASInfo"]][["varComp"]]), 2)
expect_inherits(mtg0[["GWAResult"]], "list")
expect_equal(length(mtg0[["GWAResult"]]), 1)
expect_equal(names(mtg0[["GWAResult"]]), "Y")

## Check results.

expect_equal_to_reference(mtg0[["GWAResult"]][["Y"]], "mtg0_GWARes")
expect_equal_to_reference(result1[["Y"]], "result1_GWARes")
expect_equal_to_reference(result2[["Y"]], "result2_GWARes")
expect_equal_to_reference(result3[["Y"]], "result3_GWARes")

## Check that specifying traits by number functions correctly.
expect_equal(result1, result1a)

## Check that specifying covar by number functions correctly.
expect_equal(result1, result1b)

### Test runMultiTraitGwas with chromosome specific kinship.

mtgM0 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                           GLSMethod = "multi")

resultM1 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              GLSMethod = "multi", covar = "V1")[["GWAResult"]]
resultM2 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              covModel = "pw", GLSMethod = "multi",
                              snpCov = "M002")[["GWAResult"]]
resultM3 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              covModel = "pw", covar = "V1",
                              GLSMethod = "multi",
                              snpCov = "M002")[["GWAResult"]]

## Check output structure.

expect_inherits(mtgM0, "GWAS")
expect_equal(length(mtgM0), 5)
expect_equal(names(mtgM0),
             c("GWAResult", "signSnp", "kinship", "thr", "GWASInfo"))
expect_equal(names(mtgM0[["GWASInfo"]]),
             c("call", "MAF", "thrType", "GLSMethod", "covModel",
               "varComp", "genomicControl", "inflationFactor"))
expect_inherits(mtgM0[["kinship"]], "list")
expect_equal(length(mtgM0[["kinship"]]), 2)
expect_inherits(mtgM0[["kinship"]][[1]], "matrix")
expect_equal(length(mtgM0[["GWASInfo"]][["varComp"]]), 2)
expect_equal(length(mtgM0[["GWASInfo"]][["varComp"]][[1]]), 2)
expect_inherits(mtgM0[["GWAResult"]], "list")
expect_equal(length(mtgM0[["GWAResult"]]), 1)
expect_equal(names(mtgM0[["GWAResult"]]), "Y")
expect_equal(mtgM0[["GWASInfo"]][["GLSMethod"]] , "multi")

## Check results.

expect_equal_to_reference(mtgM0[["GWAResult"]][["Y"]], "mtgM0_GWARes")
expect_equal_to_reference(resultM1[["Y"]], "resultM1_GWARes")
expect_equal_to_reference(resultM2[["Y"]], "resultM2_GWARes")
expect_equal_to_reference(resultM3[["Y"]], "resultM3_GWARes")


### Check option thrType

## Input checks.
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fixed", LODThr = "a"),
             "LODThr should be a single numerical value greater than 0")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "small", nSnpLOD = "a"),
             "nSnpLOD should be a single numerical value greater than 0")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fdr", alpha = "a"),
             "alpha should be a single numerical value greater than 0")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fdr", rho = "a"),
             "rho should be a single numerical value between 0 and 1")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fdr", pThr = "a"),
             "pThr should be a single numerical value between 0 and 1")

mtg_fixed <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fixed", LODThr = 1)
expect_equal(nrow(mtg_fixed$signSnp$Y), 14)

mtg_small <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "small", nSnpLOD = 5)
expect_equal(nrow(mtg_small$signSnp$Y), 10)

mtg_fdr <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               thrType = "fdr")
expect_equal(nrow(mtg_fdr$signSnp$Y), 6)


### Check option fitVarComp

expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               fitVarComp = FALSE),
             "Vg should be a matrix")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               fitVarComp = FALSE, Vg = Vg),
             "Ve should be a matrix")

Vg0 <- Vg
Ve0 <- Ve
colnames(Vg0)[1] <- colnames(Ve0)[1] <- "t4"

expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               fitVarComp = FALSE, Vg = Vg0, Ve = Ve),
             "Column names and rownames of Vg should be identical and")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               fitVarComp = FALSE, Vg = Vg, Ve = Ve0),
             "Column names and rownames of Ve should be identical and")

mtg_VgVe <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              fitVarComp = FALSE, Vg = Vg, Ve = Ve)
expect_equivalent(mtg_VgVe$GWASInfo$varComp$Vg, Vg[1:2, 1:2])
expect_equivalent(mtg_VgVe$GWASInfo$varComp$Ve, Ve[1:2, 1:2])


### Check option estCom

mtg_estCom <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                                estCom = TRUE)

expect_equal(mtg0$GWAResult$Y, mtg_estCom$GWAResult$Y[, 1:9])

expect_equal(colnames(mtg_estCom$GWAResult$Y)[10:13],
             c("pValCom", "effsCom", "effsComSe", "pValQtlE"))

expect_equal_to_reference(mtg_estCom$GWAResult$Y, "result_estCom")


### Check options mG and mE

expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               covModel = "fa", mG = "a"),
             "mG should be a single numerical value between 1 and 1")
expect_error(runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                               covModel = "fa", mE = "a"),
             "mE should be a single numerical value between 1 and 1")

mtg_mGmE <- runMultiTraitGwas(gData = gDataTest, covModel = "fa",
                              mG = 2, mE = 2, maxIter = 100)

expect_equal_to_reference(mtg_mGmE$GWAResult$Y, "result_mGmE")



