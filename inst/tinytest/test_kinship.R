## Test multi allele kinship function.

set.seed(1234)
map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                  pos = c(1, 4, 17, 2, 3, 5))
posCor <- c(3, 8, 13, 1, 1.5, 2)
X <- matrix(runif(n = 180), nrow = 60)
X <- X / rowSums(X)
dim(X) <- c(10, 6, 3)
K0 <- statgenQTLxT:::multiAllKin(x = X, posCor = posCor)

expect_true(inherits(K0, "matrix"))
expect_equal(dim(K0), c(10, 10))
expect_equal(diag(K0), rep(x = 1, times = 10))
expect_equal(K0[, 1],
             c(1, 0.318199140315202, 0.353188959968734, 0.356389520929043,
               0.365703566465929, 0.344204968497873, 0.363187385587088,
               0.353630750320603, 0.303310956779327, 0.383196299221719))
K1 <- statgenQTLxT:::multiAllKin(x = X, posCor = posCor, denom = 1)
expect_equal(K1 / 28.5, K0)

## Test within kinship function.

K2 <- kinship(X = X, map = map, method = "multiAllKin")
expect_equal(K2, K0)
K3 <- kinship(X = X, map = map, method = "multiAllKin", denominator = 1)
expect_equal(K3, K1)


