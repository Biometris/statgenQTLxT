load(file = "testdata.rda")

gDataTest <- statgenGWAS::createGData(geno = X, map = map, pheno = Y, covar = W)

### Test check functions.

## Test chkGData.
expect_error(statgenQTLxT:::chkGData(), "valid gData object")
expect_error(statgenQTLxT:::chkGData(X), "valid gData object")
expect_silent(statgenQTLxT:::chkGData(gDataTest))

# Copy and remove markers to create an invalid gData object.
gDataTest2 <- gDataTest
gDataTest2$markers <- NULL
expect_error(statgenQTLxT:::chkGData(gDataTest2), "should contain markers")

# Map is still present so if only map is needed check should succeed.
expect_silent(statgenQTLxT:::chkGData(gDataTest2, comps = "map"))


## Test chkTraits.

# Input should be character or numeric.
expect_error(statgenQTLxT:::chkTraits(traits = TRUE, trials = "Y",
                                      gData = gDataTest, multi = TRUE),
             "should be a numeric or character vector")

# When multi = FALSE only single numeric/character values are accepted.
expect_error(statgenQTLxT:::chkTraits(traits = 1:2, trials = "Y",
                                      gData = gDataTest, multi = FALSE),
             "should be a single numeric or character value")

# Check character input.
expect_silent(statgenQTLxT:::chkTraits(traits = "t1", trials = "Y",
                                       gData = gDataTest, multi = FALSE))
expect_silent(statgenQTLxT:::chkTraits(traits = c("t1", "t2"), trials = "Y",
                                       gData = gDataTest, multi = TRUE))
expect_error(statgenQTLxT:::chkTraits(traits = "t6", trials = "Y",
                                      gData = gDataTest, multi = FALSE),
             "For Y not all traits")
expect_error(statgenQTLxT:::chkTraits(traits = c("t1", "t6"),
                                      trials = "Y",
                                      gData = gDataTest, multi = TRUE),
             "For Y not all traits")

# Check numeric input.
expect_silent(statgenQTLxT:::chkTraits(traits = 2, trials = 1, gData = gDataTest,
                                       multi = FALSE))
expect_silent(statgenQTLxT:::chkTraits(traits = c(2, 3), trials = 1,
                                       gData = gDataTest, multi = TRUE))
expect_error(statgenQTLxT:::chkTraits(traits = 1, trials = 1, gData = gDataTest,
                                      multi = FALSE),
             "For 1 not all traits")
expect_error(statgenQTLxT:::chkTraits(traits = c(2, 7), trials = 1,
                                      gData = gDataTest, multi = TRUE),
             "For 1 not all traits")


## Test chkTrials.

# Input should be character or numeric.
expect_error(statgenGWAS:::chkTrials(trials = TRUE, gData = gDataTest),
             "should be a numeric or character vector")

# Check character input.
expect_silent(statgenGWAS:::chkTrials(trials = "Y", gData = gDataTest))
expect_error(statgenGWAS:::chkTrials(trials = "ph3", gData = gDataTest),
             "should be in pheno")

# Check numeric input.
expect_silent(statgenGWAS:::chkTrials(trials = 1, gData = gDataTest))
expect_error(statgenGWAS:::chkTrials(trials = 3, gData = gDataTest),
             "should be in pheno")


## Test chkNum

expect_silent(statgenQTLxT:::chkNum(1, min = 0, max = 2))
expect_error(statgenQTLxT:::chkNum("a"), "single numerical value")
expect_error(statgenQTLxT:::chkNum(c(1, 2)), "single numerical value")
expect_error(statgenQTLxT:::chkNum(1, min = NULL, max = 0), "smaller than 0")
expect_error(statgenQTLxT:::chkNum(1, min = 2, max = NULL), "greater than 2")
expect_error(statgenQTLxT:::chkNum(1, min = 2, max = 3), "between 2 and 3")

## Test chkMarkers

# Input should be numerical.
testMrk <- matrix(letters[1:4], nrow = 2)
expect_error(statgenQTLxT:::chkMarkers(markers = testMrk),
             "should be a numerical matrix")

# No missing values allowed.
testMrk2 <- matrix(c(1:3, NA), nrow = 2)
expect_error(statgenQTLxT:::chkMarkers(markers = testMrk2),
             "markers contains missing values")

expect_error(statgenQTLxT:::chkMarkers(markers = testMrk, dim = 3),
             "markers should be a three-dimensional array")

expect_silent(statgenQTLxT:::chkMarkers(markers = gDataTest$markers))

## Test chkCovar

# Input should be character or numeric.
expect_error(statgenQTLxT:::chkCovar(covar = TRUE, gData = gDataTest),
             "should be a numeric or character vector")

# Covariates should be in covar in gData.
expect_error(statgenQTLxT:::chkCovar(covar = "V3", gData = gDataTest),
             "covar should be columns in covar in gData")

expect_silent(statgenQTLxT:::chkCovar(covar = "V1", gData = gDataTest))

## Test chkSnpCov

# SNP Covariates should be in markers in gData.
expect_error(statgenQTLxT:::chkSnpCov(snpCov = "SNP1", gData = gDataTest),
             "All snpCov should be in markers")

expect_silent(statgenQTLxT:::chkSnpCov(snpCov = "M001", gData = gDataTest))

## Test chkKin

# When GLSMethod = "single" kin should be a matrix.
expect_error(statgenQTLxT:::chkKin(kin = 1, gData = gDataTest,
                                   GLSMethod = "single"),
             "kin should be a matrix")
expect_silent(statgenQTLxT:::chkKin(kin = K, gData = gDataTest,
                                    GLSMethod = "single"))

# When GLSMethod = "multi" kin should be a list of matrices.
expect_error(statgenQTLxT:::chkKin(kin = K, gData = gDataTest,
                                   GLSMethod = "multi"),
             "kin should be a list of matrices of length equal to the number")
expect_silent(statgenQTLxT:::chkKin(kin = list(K, K), gData = gDataTest,
                                    GLSMethod = "multi"))

