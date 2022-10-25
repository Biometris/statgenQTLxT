load(file = "testdata.rda")

gDataTest <- statgenGWAS::createGData(geno = X, map = map, pheno = Y, covar = W)

### Test runMultiTraitGwas with single kinship.

mtg0 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"))

result1 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             covar = "V1")[["GWAResult"]]
result1a <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              covar = 1)[["GWAResult"]]
result2 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             snpCov = "M002")[["GWAResult"]]
result3 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                             covar = "V1", snpCov = "M002")[["GWAResult"]]

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

expect_equal_to_reference(mtg0[["GWAResult"]][["Y"]],
                          "mtg0_GWARes")
expect_equal_to_reference(result1[["Y"]],
                          "result1_GWARes")
expect_equal_to_reference(result2[["Y"]],
                          "result2_GWARes")
expect_equal_to_reference(result3[["Y"]],
                          "result3_GWARes")

## Check that specifying covar by number functions correctly.
expect_equal(result1, result1a)

### Test runMultiTraitGwas with chromosome specific kinship.

mtgM0 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                           GLSMethod = "multi")
resultM1 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              GLSMethod = "multi", covar = "V1")[["GWAResult"]]
resultM2 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              GLSMethod = "multi",
                              snpCov = "M002")[["GWAResult"]]
resultM3 <- runMultiTraitGwas(gData = gDataTest, traits = c("t1", "t2"),
                              covar = "V1", GLSMethod = "multi",
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

expect_equal_to_reference(mtgM0[["GWAResult"]][["Y"]],
                          "mtgM0_GWARes")
expect_equal_to_reference(resultM1[["Y"]],
                          "resultM1_GWARes")
expect_equal_to_reference(resultM2[["Y"]],
                          "resultM2_GWARes")
expect_equal_to_reference(resultM3[["Y"]],
                          "resultM3_GWARes")
