load(file = "testdata.rda")

## Test estimation of effects.

effs0 <- statgenQTLxT:::estEffsCPP(y = Y, w = W[, 1, drop = FALSE], x = X,
                                   vg = Vg, ve = Ve, k = K)
effs1 <- statgenQTLxT:::estEffsCPP(y = Y, w = W[, 1, drop = FALSE], x = X,
                                   vg = Vg, ve = Ve, k = K, estCom = TRUE)

# Check structure.

expect_true(inherits(effs0, "list"))
expect_equal(length(effs0), 7)
expect_equal(names(effs0), c("effs", "effsSe", "pVals", "effsCom", "effsComSe",
                             "pValsCom", "pValsQtlE"))
expect_equal(length(effs1), 7)
expect_equal(names(effs1), c("effs", "effsSe", "pVals", "effsCom", "effsComSe",
                             "pValsCom", "pValsQtlE"))
expect_true(inherits(effs1[["effs"]], "matrix"))
expect_true(inherits(effs1[["effsSe"]], "matrix"))
expect_equal(dim(effs1[["effs"]]), c(3, 100))
expect_equal(dim(effs1[["effsSe"]]), c(3, 100))
expect_true(inherits(effs1[["pVals"]], "numeric"))
expect_true(inherits(effs1[["effsCom"]], "numeric"))
expect_true(inherits(effs1[["effsComSe"]], "numeric"))
expect_true(inherits(effs1[["pValsCom"]], "numeric"))
expect_true(inherits(effs1[["pValsQtlE"]], "numeric"))
expect_equal(effs0[1:3], effs1[1:3])

# Check results.

expect_equal(effs1[["effs"]][1:10],
             c(0.241140292035063, 0.566426284782653, 0.200563422195687,
               0.5072980814268, 0.895706387676149, 0.329167063086988,
               0.774473488319093, 1.0813170308548, 0.450996816155276,
               0.209075775078276))
expect_equal(effs1[["effsSe"]][1:10],
             c(0.0710897574205907, 0.0995815912538212, 0.0742516062831747,
               0.0837727845895638, 0.117862960438281, 0.087082081326973,
               0.0841365173082272, 0.117566751302421, 0.0880166839367443,
               0.0699792657350049))
expect_equal(effs1[["pVals"]][1:10],
             c(0.000206528512132969, 8.49066061435212e-09, 3.14068033447528e-16,
               0.00212406748361569, 1.99005437736383e-17, 3.02515119624494e-11,
               1.16179422503003e-11, 0.026689383989729, 6.24226826016099e-12,
               3.28745335819008e-05))
expect_equal(effs1[["effsCom"]][1:10],
             c(0.227613830841915, 0.428450598330427, 0.632750988225132,
               0.196821961815739, 0.629030152176383, 0.490584048466665,
               0.593511964658955, 0.206981748536258, 0.50257289427201,
               0.343158290899609))
expect_equal(effs1[["effsComSe"]][1:10],
             c(0.0596125235420024, 0.0701238866736423, 0.0705835313652413,
               0.0587284964549322, 0.0720210053506788, 0.0744071879910642,
               0.0707325319018773, 0.0615207439135682, 0.063712562738162,
               0.0644950374681231))
expect_equal(effs1[["pValsCom"]][1:10],
             c(0.00365640202193862, 2.54510801083346e-06, 1.62733991401514e-12,
               0.0108710047451377,6.60464395970753e-12, 3.52853858497932e-07,
               4.8638544160601e-11, 0.0105629736163054, 7.56949593261243e-10,
               4.54284887693718e-05))
expect_equal(effs1[["pValsQtlE"]][1:10],
             c(0.00378841700698975, 0.000102166982623902, 2.85810581051649e-06,
               0.0166611555765428, 4.47968127261315e-08, 2.16256371137555e-06,
               0.00445316399674635, 0.261196665406218, 0.000164374892733174,
               0.0326528041804481))

# Check options returnSe.

effs0a <- statgenQTLxT:::estEffsCPP(y = Y, w = W[, 1, drop = FALSE], x = X,
                                    vg = Vg, ve = Ve, k = K,
                                    returnSe = FALSE)
effs1a <- statgenQTLxT:::estEffsCPP(y = Y, w = W[, 1, drop = FALSE], x = X,
                                    vg = Vg, ve = Ve, k = K,
                                    returnSe = FALSE, estCom = TRUE)
expect_true(all(effs0a[["effsSe"]] == 0))
expect_true(all(effs0a[["pVals"]] == 0))
expect_true(all(effs0a[["effsCom"]] == 0))
expect_true(all(effs1a[["effsSe"]] == 0))
expect_true(all(effs1a[["pVals"]] == 0))
expect_equal(effs0a[["effs"]], effs0[["effs"]])
expect_equal(effs1a[c("effs", "effsCom")], effs1[c("effs", "effsCom")])

# Add covariates.

effs0c <- statgenQTLxT:::estEffsCPP(y = Y, w = W, x = X, vg = Vg, ve = Ve,
                                    k = K)
effs1c <- statgenQTLxT:::estEffsCPP(y = Y, w = W, x = X, vg = Vg, ve = Ve,
                                    k = K, estCom = TRUE)
expect_equal(effs0c[c("effs", "effsSe", "pVals")],
             effs1c[c("effs", "effsSe", "pVals")])
expect_equal(effs1c[["effs"]][1:10],
             c(0.123234858949716, 0.362182223951438, 0.111445152593952,
               0.347721242501487, 0.606512198557773, 0.198085207680003,
               0.650236348332276, 0.816266552194951, 0.331812523023079,
               0.0987141903510034))
expect_equal(effs1c[["effsSe"]][1:10],
             c(0.0735388920667274, 0.103500018489474, 0.0765457907358536,
               0.0894045814720542, 0.127048931988649, 0.0924329576968679,
               0.0905653140277239, 0.127885230423261, 0.0940125236536251,
               0.072091806246925))
expect_equal(effs1c[["pVals"]][1:10],
             c(0.0377937210405903, 0.000349428417481923, 5.08247369928231e-10,
               0.129440922522865, 4.42703463340778e-09, 7.75713184827169e-07,
               3.97537388934867e-06, 0.223346100447953, 1.19743079516545e-06,
               0.00182584576049935))
expect_equal(effs1c[["effsCom"]][1:10],
             c(0.00290917767144171, 0.00220357049812621, 0.0019973160219841,
               0.00280647029928628, 0.00181598981880013, 0.00188028051559034,
               0.00208772594995757, 0.00293353143124234, 0.00213019286900433,
               0.0023860192211404))
expect_equal(effs1c[["effsComSe"]][1:10],
             c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN))
expect_equal(effs1c[["pValsCom"]][1:10],
             c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
expect_equal(effs1c[["pValsQtlE"]][1:10],
             c(7.14386821939662e-17, 5.17494821660331e-19, 4.98118594739604e-25,
               2.89008041179995e-16, 4.47862278987384e-24, 8.90845453872134e-22,
               5.0010678523824e-21, 5.64913691565916e-16, 1.42016279181918e-21,
               2.86167708757127e-18))
