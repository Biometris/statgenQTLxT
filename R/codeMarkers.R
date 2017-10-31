#' Code and impute markers
#'
#' \code{codeMarkers} codes markers in a \code{gData} object and optionally performs imputation
#' of missing values as well.\cr
#' The function performs the following steps:\cr
#' \enumerate{
#' \item{remove SNPs with more missing values then \code{nMiss}.}
#' \item{code SNPs to numerical values.}
#' \item{remove SNPs with a minor allele frequency lower than \code{MAF}.}
#' \item{optionally remove duplicate SNPs.}
#' \item{optionally impute missing values.}
#' \item{repeat steps 3. and 4. if missing values are imputed.}
#' }
#'
#' @param gData an object of class \code{gData} containing at least \code{markers}.
#' @param refAll a character string indicating the reference allele used when recoding markers.\cr
#' If "minor" then the recoding is done using the minor allele as reference allele. Alternatively
#' a single character can be supplied as a reference allele for the whole set of SNPs, or a
#' character vector with a reference allele per SNP.
#' @param nMiss a numical value between 0 and 1. SNPs with a fraction of missings higher then
#' \code{nMiss} will be removed. SNPs with only missing values will always be removed.
#' @param MAF a numerical value between 0 and 1. SNPs with a Minor Allele Frequency (MAF) below
#' this value will be removed.
#' @param removeDuplicates should duplicate SNPs be removed?
#' @param keep a vector of SNPs. These SNPs will never be removed in the whole process.
#' @param impute should imputation of missing values be done?
#' @param imputeType a character string indicating what kind of imputation of values should be done.\cr
#' \itemize{
#' \item{fixed - missing values will be replaced by a given fixed value.}
#' \item{random - missing values will be replaced by a random value calculated using allele
#' frequencies per SNP.}
#' \item{beagle - missing values will be imputed using beagle software. If you use this option
#' please cite the original papers in your publication.}
#' }
#' @param fixedValue numerical value used for replacing missing values in case \code{inputType}
#' is fixed.
#'
#' @return a copy of the input \code{gData} object with markers replaced by coded and imputed markers.
#'
#' @references {B L Browning and S R Browning (2016). Genotype imputation with millions of reference
#' samples. Am J Hum Genet 98:116-126. doi:10.1016/j.ajhg.2015.11.020}
#'
#' @export
#'
#' @examples ## Create markers
#' markers <- matrix(c(
#' "AA",   "AB",   "AA",   "BB",   "BA",   "AB",   "AA",   "AA",  NA, "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "BB",   "AA",  NA, "AA",
#' "AA",   "BA",   "AB",   "BB",   "AB",   "AB",   "AA",   "BB",  NA, "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "AA",   "AA",   "AA",  NA, "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   "BB",   "BB",   "BB",  "AB", "AA",
#' "AA",   "AA",   "BB",   "BB",   "AA",   NA,     "BB",   "AA",  NA, "AA",
#' "AB",   "AB",   "BB",   "BB",   "BB",   "AA",   "BB",   "BB",  NA, "AB",
#' "AA",   "AA",    NA,    "BB",    NA,    "AA",   "AA",   "AA",  "AA", "AA",
#' "AA",    NA,     NA,    "BB",   "BB",   "BB",   "BB",   "BB",  "AA", "AA",
#' "AA",    NA,    "AA",   "BB",   "BB",   "BB",   "AA",   "AA",  NA, "AA"),
#' ncol = 10, byrow = TRUE, dimnames = list(paste0("IND", 1:10), paste0("SNP", 1:10)))
#'
#' ## create object of class 'gData'.
#' gData <- createGData(geno = markers)
#'
#' ## Code markers by minor allele, no imputation.
#' gDataCoded1 <- codeMarkers(gData = gData, impute = FALSE)
#'
#' ## Code markers by reference alleles, impute missings by fixed value.
#' gDataCoded2 <- codeMarkers(gData = gData, refAll = rep(x = c("A", "B"), times =  5),
#'                            impute = TRUE, imputeType = "fixed", fixedValue = 1)
#'
#' ## Code markers by minor allele, impute by random value.
#' gDataCoded3 <- codeMarkers(gData = gData, impute = TRUE, imputeType = "random")

codeMarkers <- function(gData,
                        refAll = "minor",
                        nMiss = 1,
                        MAF = NULL,
                        removeDuplicates = TRUE,
                        keep = NULL,
                        impute = TRUE,
                        imputeType = "random",
                        fixedValue = NULL) {
  ## Checks.
  if(missing(gData) || !is.gData(gData) || is.null(gData$markers)) {
    stop("gData should be a valid gData object with at least markers included.\n")
  }
  if (length(refAll) > 1 && !length(refAll) == ncol(gData$markers)) {
    stop("number of reference alleles should either be 1 or equal to the amount
         of SNPs in markers.\n")
  }
  if (is.null(nMiss) || length(nMiss) > 1 || !is.numeric(nMiss) ||
      !dplyr::between(x = nMiss, left = 0, right = 1)) {
    stop("nMiss should be a single numerical value between 0 and 1.\n")
  }
  if (!is.null(MAF) && (length(MAF) > 1 || !is.numeric(MAF) ||
                        !dplyr::between(x = MAF, left = 0, right = 1))) {
    stop("MAF should be a single numerical value between 0 and 1.\n")
  }
  if (!is.null(keep) && (!is.character(keep) || !all(keep %in% colnames(gData$markers)))) {
    stop("all items in keep should be SNPs in markers.\n")
  }
  if (impute) {
    if (length(imputeType) > 1 && !imputeType %in% c("fixed", "random", "beagle")) {
      stop("imputeType should be one of fixed, random or beagle.\n")
    }
    if (imputeType == "fixed" && is.null(fixedValue)) {
      stop("fixedValue cannot be NULL.\n")
    }
  }
  markersOrig <- as.matrix(gData$markers)
  snpKeep <- colnames(markersOrig) %in% keep
  ## Remove markers with too many missings.
  if (!is.null(nMiss)) {
    snpMiss <- colMeans(is.na(markersOrig)) <= nMiss
    markersClean <- markersOrig[, snpMiss | snpKeep]
    snpKeep <- snpKeep[which(snpMiss | snpKeep)]
    if (length(refAll) > 1) {
      refAll <- refAll[which(snpMiss | snpKeep)]
    }
  }
  ## Recode markers.
  if (!is.numeric(markersClean)) {
    if (refAll[1] == "minor") {
      ## Set first allele per SNP as reference allele.
      alleles <- lapply(X = as.data.frame(markersClean), FUN = levels)
      firstAlls <- sapply(X = alleles, FUN = "[[", 1)
      ## Split either on puntuation or on single letters.
      refAlls <- sapply(X = strsplit(x = firstAlls,
                                     split = ifelse(length(grep(pattern = "[[:punct:]]",
                                                                x = firstAlls)) == 0,
                                                    "", "[[:punct:]]")),
                        FUN = "[[", 1)
      ## Expand refAll to match dimensions markersClean.
      refAlls <- rep(x = refAlls, each = nrow(markersClean))
      ## Get maximum number of alleles.
      maxAll <- unname(nchar(gsub(pattern = "[[:punct:]]",
                                  replacement = "",
                                  x = alleles[[1]][1])) %/% nchar(refAlls[1]))
    } else if (length(refAll) > 1) {
      ## Expand refAll to match dimensions markersClean.
      refAlls <- rep(refAll, each = nrow(markersClean))
      ## Get maximum number of alleles.
      maxAll <- nchar(markersClean[ , 1][!is.na(markersClean[, 1])][1]) %/% nchar(refAll[1])
    } else {
      refAlls <- refAll
      ## Get maximum number of alleles.
      maxAll <- nchar(markersClean[ , 1][!is.na(markersClean[, 1])][1]) %/% nchar(refAll)
    }
    ## Recode markers by counting the number of times the reference marker is in the alleles.
    ## nchar keeps NA so no need to work around NAs.
    markersRecoded <- matrix(mapply(FUN = function(x, y) {
      (nchar(x) - nchar(gsub(pattern = y, replacement = "", x = x))) / nchar(y)
    }, markersClean, refAlls, USE.NAMES = FALSE),
    nrow = nrow(markersClean),
    dimnames = dimnames(markersClean))
    if (refAll[1] == "minor") {
      ## Correct for position of minor allele.
      markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2] <-
        maxAll - markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2]
    }
  }
  ## Remove markers with low MAF.
  if (!is.null(MAF)) {
    snpMAF <- dplyr::between(x = colMeans(markersRecoded, na.rm = TRUE),
                             left = maxAll * MAF, right = maxAll * (1 - MAF))
    markersRecoded <- markersRecoded[, snpMAF | snpKeep]
    snpKeep <- snpKeep[which(snpMAF | snpKeep)]
  }
  ## Remove duplicated markers.
  if (removeDuplicates) {
    if (anyDuplicated(markersRecoded, MARGIN = 2)) {
      ## Only using duplicated would always remove the first occurence.
      ## Using sample to make it random, always putting keep SNPs first.
      randOrder <- c(c(1:ncol(markersRecoded))[snpKeep],
                     sample(x = c(1:ncol(markersRecoded))[!snpKeep]))
      dubs <- duplicated(markersRecoded[, randOrder], MARGIN = 2)[order(randOrder)]
      markersRecoded <- markersRecoded[, !dubs | snpKeep]
      snpKeep <- snpKeep[which(!dubs | snpKeep)]
    }
  }
  ## Impute missing values.
  if (impute) {
    if (imputeType == "fixed") {
      ## Replace missing values by fixed value.
      markersRecoded[is.na(markersRecoded)] <- fixedValue
    } else if (imputeType == "random") {
      ## Replace missing values by random value based on probabilities per SNP.
      snpNA <- apply(X = markersRecoded, MARGIN = 2, FUN = anyNA)
      markersRecoded[, snpNA] <- apply(X = markersRecoded[, snpNA], MARGIN = 2,
                                       FUN = function(x) {
                                         p <- mean(x, na.rm = TRUE) / maxAll
                                         x[is.na(x)] <- sample(x = 0:maxAll,
                                                               size = sum(is.na(x)),
                                                               replace = TRUE,
                                                               prob = choose(maxAll, 0:maxAll) *
                                                                 (1 - p) ^ (maxAll:0) *
                                                                 p ^ (0:maxAll))
                                         return(x)
                                       })
    } else if (imputeType == "beagle") {
      ## Imputation of missing values using beagle software.
      if (!dir.exists("beagle")) {
        dir.create("beagle")
      }
      ## Set prefix for distinguishing file names.
      prefix <- format(Sys.time(), "%y%m%d%H%M%OS3")
      map <- gData$map[colnames(markersRecoded), ]
      ## Convert map to format suitable for beagle input.
      mapBeagle <- data.frame(map$chr,
                              rownames(map),
                              map$pos,
                              map$pos)
      if (!is.integer(mapBeagle[, 1])) {
        mapBeagle[, 1] <- as.integer(as.factor(mapBeagle[, 1]))
      }
      mapBeagle[, 1] <- paste0("chr", mapBeagle[, 1])
      ## Write map to .map file
      write.table(mapBeagle, file = paste0("beagle/run",
                                           prefix, ".map"), col.names = FALSE,
                  row.names = FALSE, quote = FALSE, na = ".",
                  sep = "\t")
      ## Convert markers to format suitable for beagle input.
      all00 <- "0/0"
      all01 <- "0/1"
      all11 <- "1/1"
      all10 <- "1/0"
      markersBeagle <- as.data.frame(t(markersRecoded), stringsAsFactors = FALSE)
      markersBeagle[markersBeagle == 0] <- all00
      markersBeagle[markersBeagle == 1] <- all01
      markersBeagle[markersBeagle == 2] <- all11
      markersBeagle[markersBeagle == -1] <- all10
      ## Write markers to vcf file.
      vcfBeagle <- cbind(data.frame(CHROM = mapBeagle[, 1], POS = mapBeagle[, 3],
                                    ID = rownames(map), REF = "A", ALT = "G", QUAL = ".",
                                    FILTER = "PASS", INFO = ".", FORMAT = "GT",
                                    stringsAsFactors = FALSE),
                         markersBeagle)
      vcfFile = paste0("beagle/run", prefix, "input.vcf")
      cat(file = vcfFile,
          "##fileformat=VCFv4.1",
          "\n##filedate=", Sys.Date(),
          "\n##source=\"codeMarkers of genStatPipeline\"",
          "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          "\n#")
      cat(file = vcfFile, paste(colnames(vcfBeagle), collapse = "\t"), "\n",
          append = TRUE)
      write.table(vcfBeagle, file = vcfFile, quote = FALSE, col.names = FALSE,
                  row.names = FALSE, append = TRUE, sep = "\t",
                  na = paste(".", ".", sep = "/"))
      ## Run beagle with default settings.
      system(paste0("java -Xmx3000m -jar ",
                    shQuote(paste0(sort(path.package()[grep("genStatPipeline",
                                                            path.package())])[1],
                                   "/java/beagle.jar")), " gtgl=beagle/run",
                    prefix, "input.vcf out=beagle/run",
                    prefix, "out gprobs=true nthreads=", 1,
                    " map=beagle/run", prefix, ".map "), intern = TRUE)
      ## Read beagle output.
      beagleOut <- read.table(gzfile(paste0("beagle/run", prefix, "out.vcf.gz")),
                              stringsAsFactors = FALSE)
      ## Convert beagle output to format suitable for gData.
      markersRecoded <- t(beagleOut[, 10:ncol(beagleOut)])
      markersRecoded[substr(markersRecoded, 1, 3) == all00] <- 0
      markersRecoded[substr(markersRecoded, 1, 3) == all01] <- 1
      markersRecoded[substr(markersRecoded, 1, 3) == all11] <- 2
      markersRecoded[substr(markersRecoded, 1, 3) == all10] <- NA
      markersRecoded <- apply(X = markersRecoded, MARGIN = 2, FUN = as.numeric)
      rownames(markersRecoded) <- colnames(markersBeagle)
      colnames(markersRecoded) <- beagleOut[, 3]
    }
    if (refAll[1] == "minor") {
      ## Correct for position of minor allele after imputation.
      markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2] <-
        maxAll - markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > maxAll / 2]
    }
    ## Remove markers with low MAF after imputation.
    if (!is.null(MAF)) {
      snpMAF <- dplyr::between(colMeans(markersRecoded, na.rm = TRUE),
                               maxAll * MAF, maxAll * (1 - MAF))
      markersRecoded <- markersRecoded[, snpMAF | snpKeep]
      snpKeep <- snpKeep[which(snpMAF | snpKeep)]
    }
    ## Remove duplicated markers after imputation.
    if (anyDuplicated(markersRecoded)) {
      randOrder <- c(c(1:ncol(markersRecoded))[snpKeep],
                     sample(x = c(1:ncol(markersRecoded))[!snpKeep]))
      dubs <- duplicated(markersRecoded[, randOrder], MARGIN = 2)[order(randOrder)]
      markersRecoded <- markersRecoded[, !dubs]
      snpKeep <- snpKeep[which(!dubs | snpKeep)]
    }
  }
  ## Return gData object with recoded and imputed markers.
  ## SuppressWarnnings needed since original geno will be overwritten.
  return(suppressWarnings(createGData(gData = gData, geno = markersRecoded)))
}





