## Create data for testing functions.
## Restricted and anonymized version of DROPs data.
#load(system.file("extdata", "testData.RData", package = "statgenPipeline"))

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

X <- t(t(sample(x = 2, size = 100, replace = TRUE)))

rownames(Y) <- rownames(X) <- rownames(K) <- colnames(K) <-
  paste0("G", formatC(1:100, width = 3, flag = "0"))
colnames(Y) <- paste0("t", 1:3)
colnames(X) <- "V1"


# ## Create a dataset for unit testing.
# set.seed(1234)
# y <- matrix(rnorm(50, mean = 10, sd = 2), nrow = 10)
# X <- matrix(sample(x = c(0, 1), size = 30, replace = TRUE), nrow = 10)
# Sigma <- matrix(runif(n = 100), nrow = 10)
# ## Assure sigma is symmetric.
# Sigma <- tcrossprod(Sigma)
# covs <- matrix(runif(n = 20, max = 100), nrow = 10)
# pheno <- data.frame(genotype = paste0("G", 1:10), y)
# ## Add random missing values to pheno2.
# pheno2 <- pheno
# for (i in 2:6) {
#   pheno2[sample(x = 1:10, size = 2), i] <- NA
# }
# map <- data.frame(chr = c(1, 1, 2), pos = 1:3)
# rownames(X) <- rownames(y) <- rownames(Sigma) <- colnames(Sigma) <-
#   rownames(covs) <- paste0("G", 1:10)
# colnames(X) <- rownames(map) <- paste0("M", 1:3)
# colnames(y) <- paste0("t", 1:5)
# ## Create gData object.
# gDataTest <- createGData(map = map, geno = X, kin = Sigma,
#                          pheno = list(ph1 = pheno, ph2 = pheno2),
                         # covar = as.data.frame(covs))
## Export to package
save(#gDataTest, X, y, Sigma, map, covs,
     Y, K, X,
     file = "inst/tinytest/testdata.rda")
