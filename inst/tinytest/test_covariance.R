load(file = "testdata.rda")

### Test computation of covariance matrices.

## Test unstructured covariance.

covUnst <- statgenQTLxT:::covUnstr(Y = Y, K = K)

# Check structure.

expect_true(inherits(covUnst, "list"))
expect_equal(length(covUnst), 2)
expect_equal(names(covUnst), c("Vg", "Ve"))
expect_true(inherits(covUnst[["Vg"]], "matrix"))
expect_true(inherits(covUnst[["Ve"]], "matrix"))

# Check results.

expect_equal(as.numeric(covUnst[["Vg"]]),
             c(0.192371054661003, 0.207143343370603, 0.0695190613823293,
               0.207143343370603, 0.562742860825983, 0.150656062326219,
               0.0695190613823293, 0.150656062326219, 0.137965708778342))
expect_equal(as.numeric(covUnst[["Ve"]]),
             c(0.218260911481604, 0.0181226134186007, 0.0823614626198486,
               0.0181226134186007, 0.163824982416043, 0.151835604456582,
               0.0823614626198486, 0.151835604456582, 0.350573161696034))

# Add covariates.

covUnstCov <- statgenQTLxT:::covUnstr(Y = Y, K = K, X = W)
expect_equal(as.numeric(covUnstCov[["Vg"]]),
             c(0.196165320462335, 0.213742429599169, 0.0719881355611577,
               0.213742429599169, 0.503904252678917, 0.167546528651886,
               0.0719881355611577, 0.167546528651886, 0.140233886248504))
expect_equal(as.numeric(covUnstCov[["Ve"]]),
             c(0.221931815403061, 0.010764425209951, 0.0835567504378285,
               0.010764425209951, 0.249355180564012, 0.129620754183919,
               0.0835567504378285, 0.129620754183919, 0.359202021730626))

# Check option fixDiag.

expect_warning(statgenQTLxT:::covUnstr(Y = Y, K = K, fixDiag = TRUE),
               "not implemented yet")

# Check option veDiag.

covUnstVeD <- statgenQTLxT:::covUnstr(Y = Y, K = K, VeDiag = TRUE)
expect_equal(as.numeric(covUnstVeD[["Vg"]]),
             c(0.218449152375622, 0.217596220872246, 0.117432384935249,
               0.217596220872246, 0.465832867214789, 0.23922669614095,
               0.117432384935249, 0.23922669614095, 0.167245047121722))
expect_equal(as.numeric(covUnstVeD[["Ve"]]),
             c(0.176667883005282, 0, 0, 0, 0.306367553950968, 0, 0, 0,
               0.307589174602025))

## Test pariwise covariance.

covPw <- statgenQTLxT:::covPW(Y = Y, K = K)

# Check structure.

expect_true(inherits(covPw, "list"))
expect_equal(length(covPw), 2)
expect_equal(names(covPw), c("Vg", "Ve"))
expect_true(inherits(covPw[["Vg"]], "matrix"))
expect_true(inherits(covPw[["Ve"]], "matrix"))

# Check results

expect_equal(as.numeric(covPw[["Vg"]]),
             c(0.187447851984919, 0.194760087881847, 0.0820347302529487,
               0.194760087881847, 0.468125200148636, 0.132403935233835,
               0.0820347302529487, 0.132403935233835, 0.145460135908064))
expect_equal(as.numeric(covPw[["Ve"]]),
             c(0.221333965843665, 0.029982530895221, 0.0591118615515787,
               0.029982530895221, 0.294276821898124, 0.1733456475758,
               0.0591118615515787, 0.1733456475758, 0.336535659103855))

# Add covariates.

covPwCov <- statgenQTLxT:::covPW(Y = Y, K = K, X = W)
expect_equal(as.numeric(covPwCov[["Vg"]]),
             c(0.186173926875856, 0.20032684419104, 0.0795057155921952,
               0.20032684419104, 0.457805148681881, 0.127690509266447,
               0.0795057155921952, 0.127690509266447, 0.149062907043252))
expect_equal(as.numeric(covPwCov[["Ve"]]),
             c(0.232945276973869, 0.0280829970920944, 0.0675505473523137,
               0.0280829970920944, 0.316970236624609, 0.193430765106595,
               0.0675505473523137, 0.193430765106595, 0.342038967381623))

# Check option veDiag.

expect_warning(statgenQTLxT:::covPW(Y = Y, K = K, fixDiag = TRUE),
               "not implemented yet")

# Check option corMat.

covPwCor <- statgenQTLxT:::covPW(Y = Y, K = K, corMat = TRUE)
expect_equal(as.numeric(covPwCor[["Vg"]]),
             c(1, -0.339741087256352, -0.632112286967153, -0.339741087256352, 1,
               -0.514031814619155, -0.632112286967153, -0.514031814619155, 1))
expect_equal(as.numeric(covPwCor[["Ve"]]),
             c(1, -0.634325377374997, -0.572869325913341, -0.634325377374997, 1,
               -0.270255571010964, -0.572869325913341, -0.270255571010964, 1))
