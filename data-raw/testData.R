## Create data for testing functions.
set.seed(1234)

trials <- c("Bol12R", "Bol12W", "Cam12R")
genotypes <- c("11430", "A3", "A310", "A347", "A374", "A375", "A554", "AS5707",
               "B100", "B104", "B105", "B106", "B107", "B108", "B109", "B110",
               "B113", "B37", "B73", "B84", "B89", "B97", "B98", "C103", "CO109",
               "CR1Ht", "D09", "DE811", "DK4676A", "DK78010", "DK78371A", "DKFAPW",
               "DKFBHJ", "DKIBO2", "DKMBST", "EA1027", "EA1163", "EA3076", "EC136",
               "EC140", "EC151", "EC169", "EC175", "EC232", "EC242C", "EC334",
               "EP10", "EP2008-18", "EP2008-22", "EP29", "EP51", "EP52", "EP55",
               "EP67", "EP72", "EP77", "EZ11A", "EZ18", "EZ31", "EZ35", "EZ36",
               "EZ37", "EZ38", "EZ40", "EZ42", "EZ47", "EZ48", "F04401", "F04402",
               "F04701", "F04702", "F05101", "F05404", "F1808", "F1890", "F218",
               "F252", "F353", "F354", "F608", "F618", "F7001", "F7019", "F7025",
               "F7028", "F7057", "F7058", "F7081", "F7082", "F712", "F748",
               "F752", "F838", "F874", "F888", "F894", "F908", "F912", "F918",
               "F922")

dropsPheno <- statgenGWAS::dropsPheno
dropsPheno <- dropsPheno[dropsPheno[["Experiment"]] %in% trials &
                           dropsPheno[["Variety_ID"]] %in% genotypes,
                         c("Experiment", "Variety_ID", "grain.yield")]
dropsPheno <- droplevels(dropsPheno)
Y <- tapply(X = dropsPheno[["grain.yield"]],
            INDEX = dropsPheno[c("Variety_ID", "Experiment")], FUN = I)

dropsMarkers <- statgenGWAS::dropsMarkers
dropsMarkers <- as.matrix(dropsMarkers[dropsMarkers[["Ind"]] %in% genotypes, -1])
K <- kinship(dropsMarkers)

X <- dropsMarkers[, 1:100]
W <- matrix(sample(x = 2, size = 200, replace = TRUE), ncol = 2)

rownames(Y) <- rownames(X) <- rownames(K) <- colnames(K) <- rownames(W) <-
  paste0("G", formatC(1:100, width = 3, flag = "0"))
colnames(Y) <- paste0("t", 1:3)
colnames(X) <- paste0("M", formatC(1:100, width = 3, flag = "0"))
colnames(W) <- c("V1", "V2")

covUnst <- statgenPipeline:::covUnstr(Y = Y, K = K)

Vg <- covUnst[["Vg"]]
Ve <- covUnst[["Ve"]]

## Export to package.
save(Y, K, X, W, Vg, Ve,
     file = "inst/tinytest/testdata.rda")
