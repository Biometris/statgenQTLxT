load(file = "testdata.rda")

### Test util functions.

## Test dfbind.
df1 <- data.frame(a = 1:2, b = 1:2)
df2 <- data.frame(a = 1:2, c = 1:2)
df3 <- data.frame(c = 1:2, d = 1:2)

# Test that output columns are correct.
expect_equal(colnames(statgenPipeline:::dfBind(list(df1, df1))), c("a", "b"))
expect_equal(colnames(statgenPipeline:::dfBind(list(df1, df2))),
             c("a", "b", "c"))
expect_equal(colnames(statgenPipeline:::dfBind(list(df1, df3))),
             c("a", "b", "c", "d"))
expect_equal(colnames(statgenPipeline:::dfBind(list(df1, df2, df3))),
             c("a", "b", "c", "d"))

# Test that NAs are inserted for missing columns.
expect_equivalent(unlist(statgenPipeline:::dfBind(list(df1, df2))),
                  c(1, 2, 1, 2, 1, 2, NA, NA, NA, NA, 1, 2))
expect_equivalent(unlist(statgenPipeline:::dfBind(list(df1, df2, df1))),
                  c(1, 2, 1, 2, 1,2, 1, 2, NA, NA, 1, 2, NA, NA, 1, 2, NA, NA))

# Test that empty data.frames are ignored when binding.

expect_equal(statgenPipeline:::dfBind(list(data.frame(), df1)), df1)
expect_equal(statgenPipeline:::dfBind(list(df1, data.frame())), df1)
expect_equal(statgenPipeline:::dfBind(list(data.frame())), data.frame())


## Test matrixRoot
M1 <- matrix(c(1:2, 2:1), nrow = 2)
expect_equal(as.numeric(statgenPipeline:::matrixRoot(crossprod(M1))),
             c(2, 1, 1, 2))

## Test reduceKinship

K1 <- statgenPipeline:::reduceKinship(K = K, nPca = 2)
expect_true(inherits(K1, "matrix"))
expect_equal(dim(K1), dim(K))
expect_equal(K1[1:10],
             c(0.275958592649926, 0.0755268623583533, 0.0871087757501913,
               0.059879262726509, 0.0617597336440798, 0.0316564097963844,
               0.0913098771271908, -0.0120651280191346, 0.0496193372132824,
               -0.190509001746542))

## Test nearestPD

set.seed(1234)
M <- matrix(data = runif(9), nrow = 3)
## Assure symmetry.
M <- M + t(M)

M1 <- statgenPipeline:::nearestPD(M)
expect_true(inherits(M1, "matrix"))
expect_true(isSymmetric(M1))
expect_equal(as.numeric(M1),
             c(0.604727709553899, 1.04559967845123, 0.586415044934283,
               1.04559967845123, 1.82792523675371, 0.890017983604708,
               0.586415044934283, 0.890017983604709, 1.33494200875219))

# Check parameters in nearestPD.

M2 <- statgenPipeline:::nearestPD(M, corr = TRUE)
expect_equal(as.numeric(M2),
             c(1, 0.984775432717047, 0.680679114763, 0.984775432717047, 1,
               0.7976616570699, 0.680679114763, 0.7976616570699, 1))

expect_warning(statgenPipeline:::nearestPD(M, keepDiag = TRUE),
               "did not converge in 100 iterations")
M3 <- statgenPipeline:::nearestPD(M, keepDiag = TRUE, maxIter = 120)
expect_equal(as.numeric(M3),
             c(0.227406822610646, 0.614704685973532, 0.414723676348398,
               0.614704685973532, 1.72183076711372, 0.934827342556341,
               0.414723676348398, 0.934827342556341, 1.3321675164625))

M4 <- statgenPipeline:::nearestPD(M, do2eigen = FALSE)
expect_equal(as.numeric(M4),
             c(0.604727709553899, 1.04559971368566, 0.586415058601556,
               1.04559971368566, 1.82792523675371, 0.890017984212184,
               0.586415058601556, 0.890017984212184, 1.33494200875219))

M5 <- statgenPipeline:::nearestPD(M, doDykstra = FALSE)
expect_equal(as.numeric(M5),
             c(0.604727709553899, 1.04559967845123, 0.586415044934283,
               1.04559967845123, 1.82792523675371, 0.890017983604708,
               0.586415044934283, 0.890017983604708, 1.33494200875219))

M6 <- statgenPipeline:::nearestPD(M, doSym = TRUE)
expect_equal(as.numeric(M6),
             c(0.604727709553899, 1.04559967845123, 0.586415044934283,
               1.04559967845123, 1.82792523675371, 0.890017983604708,
               0.586415044934283, 0.890017983604708, 1.33494200875219))

