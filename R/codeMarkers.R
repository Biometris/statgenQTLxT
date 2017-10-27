#' Code and impute markers
#'
#' \code{codeMarkers} codes markers in a \code{gData} object and optionally performs imputation
#' of missing values as well.\cr
#' The function performs the following steps:\cr
#' \enumerate{
#' \item{remove SNPs with more missing values then \code{nMiss}.}
#' \item{code SNPs to numerical values.}
#' \item{remove SNPs with a minor allele frequency lower than \code{MAF}.}
#' \item{remove duplicate SNPs.}
#' \item{optionally impute missing values.}
#' \item{repeat steps 3. and 4. if mssing values are imputed.}
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
#' @param fixValue numerical value used for replacing missing values in case \code{inputType}
#' is fixed.
#'
#' @return the input \code{gData} object with markers replaced by coded and imputed markers.
#'
#' @references {B L Browning and S R Browning (2016). Genotype imputation with millions of reference
#' samples. Am J Hum Genet 98:116-126. doi:10.1016/j.ajhg.2015.11.020}
#'
#' @export

codeMarkers <- function(gData,
                        refAll = "minor",
                        nMiss = 1,
                        MAF = NULL,
                        keep = NULL,
                        impute = TRUE,
                        imputeType = "fix",
                        fixValue = NULL) {
  ## Checks.
  if(missing(gData) || !is.gData(gData) || is.null(gData$markers)) {
    stop("gData should be a valid gData object with at least markers included.\n")
  }
  if (is.null(nMiss) || length(nMiss) > 1 || !is.numeric(nMiss) || nMiss < 0 || nMiss > 1)
    stop("nMiss should be a single numerical value between 0 and 1.\n")
  markersOrig <- as.matrix(gData$markers)
  snpKeep <- colnames(markersOrig) %in% keep
  ## Remove markers with too many missings.
  if (!is.null(nMiss)) {
    snpMiss <- colMeans(is.na(markersOrig)) <= nMiss
    markersClean <- markersOrig[, snpMiss | snpKeep]
    snpKeep <- snpKeep[which(snpMiss)]
    if (length(refAll) > 1) {
      refAll <- refAll[which(snpMiss)]
    }
  }
  ## Recode markers.
  if (refAll[1] == "minor") {
    ## Get alleles per SNP.
    if (is.numeric(markersClean)) {
      alleles <- apply(X = markersClean, MARGIN = 2, FUN = unique, nmax = 3)
    } else {
      alleles <- lapply(X = as.data.frame(markersClean), FUN = levels)
    }
    ## Convert markers to numeric matrix.
    markersRecoded <- data.matrix(as.data.frame(markersClean))
    ## Get position of heterozygous allele per SNP.
    hetAllelePos <- unlist(lapply(X = alleles, FUN = function(x) {
      which(substr(x, 1, 1) != substr(x, nchar(x), nchar(x)))
    }))
    ## Recode entries per SNP for SNPs with heterozygous allele on position 1.
    markersRecoded[, names(which(hetAllelePos == 1))] <-
      unlist(lapply(X = markersRecoded[, names(which(hetAllelePos == 1))],
                    FUN = function (x) {
                      2 * (x %% 2) + (x %/% 2)
                    }))
    ## Recode entries per SNP for SNPs with heterozygous allele on position 3.
    markersRecoded[, names(which(hetAllelePos == 3))] <-
      unlist(lapply(X = markersRecoded[, names(which(hetAllelePos == 3))],
                    FUN = function (x) {
                      2 * ((x + 1) %% 2) + ((x + 1) %/% 2)
                    }))
    ## Subtract 1 to get 0-based coding.
    markersRecoded <- markersRecoded - 1
    ## In case there was no heterozygous allele multiply values by 2 to get 0,2 coding
    ## instead of 0,1.
    markersRecoded[ , !colnames(markersRecoded) %in% names(hetAllelePos)] <-
      2 * markersRecoded[ , !colnames(markersRecoded) %in% names(hetAllelePos)]
    ## Correct for position of minor allele.
    markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > 1] <-
      2 - markersRecoded[, colMeans(markersRecoded, na.rm = TRUE) > 1]
  } else {
    if (length(refAll) > 1) {
      ## Expand refAll to include all non-missing values in markersClean.
      refAll <- rep(refAll, each = nrow(markersClean))
    }
    ## Recode markers by counting the number of times the reference marker is in the alleles.
    ## nchar keeps NA so no need to work around NAs.
    markersRecoded <- matrix(mapply(FUN = function(x, y) {
      (nchar(x) - nchar(gsub(y, "", x))) / nchar(y)
    }, markersClean, refAll, USE.NAMES = FALSE),
    nrow = nrow(markersClean),
    dimnames = dimnames(markersClean))
  }
  ## Remove markers with low MAF.
  if (!is.null(MAF)) {
    snpMAF <- dplyr::between(colMeans(markersRecoded, na.rm = TRUE),
                             2 * MAF, 2 * (1 - MAF))
    markersRecoded <- markersRecoded[, snpMAF | snpKeep]
    snpKeep <- snpKeep[which(snpMAF)]
  }
  ## Remove duplicated markers.
  if (anyDuplicated(markersRecoded)) {
    ## Only using duplicated would always remove the first occurence.
    ## Using sample to make it random.
    randOrder <- sample(x = ncol(markersRecoded))
    dubs <- duplicated(markersRecoded[, randOrder], MARGIN = 2)[order(randOrder)]
    markersRecoded <- markersRecoded[, !dubs]
  }
  ## Impute missing values.
  if (impute) {
    if (imputeType == "fix") {
      ## Replace missing values by fixed value.
      markersRecoded[is.na(markersRecoded)] <- fixValue
    } else if (imputeType == "random") {
      ## Replace missing values by random value based on probabilities per SNP.
      snpNA <- apply(X = markersRecoded, MARGIN = 2, FUN = anyNA)
      markersRecoded[, snpNA] <- apply(X = markersRecoded[, snpNA], MARGIN = 2, FUN = function(x) {
        p <- mean(x, na.rm = TRUE) / 2
        x[is.na(x)] <- sample(x = c(0, 1, 2), size = sum(is.na(x)),
                              replace = TRUE, prob = c((1 - p) ^ 2, 2 * p, p ^ 2))
        return(x)
      })
    } else if (imputeType == "beagle") {
      ## Imputation of missing values using beagle software.
      if (!"beagle" %in% list.files()) {
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
    ## Remove markers with low MAF after imputation.
    if (!is.null(MAF)) {
      snpMAF <- dplyr::between(colMeans(markersRecoded, na.rm = TRUE),
                               2 * MAF, 2 * (1 - MAF))
      markersRecoded <- markersRecoded[, snpMAF | snpKeep]
      snpKeep <- snpKeep[which(snpMAF)]
    }
    ## Remove duplicated markers after imputation.
    if (anyDuplicated(markersRecoded)) {
      randOrder <- sample(x = ncol(markersRecoded))
      dubs <- duplicated(markersRecoded[, randOrder], MARGIN = 2)[order(randOrder)]
      markersRecoded <- markersRecoded[, !dubs]
    }
  }
  ## Return gData object with recoded and imputed markers.
  ## SuppressWarnnings needed since original geno will be overwritten.
  return(suppressWarnings(createGData(gData = gData, geno = markersRecoded)))
}





