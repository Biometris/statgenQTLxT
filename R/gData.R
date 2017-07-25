#' Create a gData object
#'
#' Create an object of S3 class gData based on external files or dataframes containing genotypic and phenotypic
#' information
#'
#' @param geno string, specifying a csv file with the genotypic data, markers in the rows and
#' genotypes in the columns. Alternatively, an dataframe with a similar layout.
#' @param genoHeader does \code{geno} contain a header row. If \code{FALSE} genotypes are given a default
#' name.
#' @param genoRowNames does \code{geno} contain rownames in its first column. If \code{FALSE} markers
#' are given a default name.
#' @param map string, specifying a csv file with at least columns for SNP name, chromosome and position.
#' The positions should be in base-pair or morgan. They should not be cumulative over the chromosomes. Other
#' columns are ignored. Alternatively, an dataframe with a similar layout.
#' @param mapSnpName the column corresponding to the SNP name in \code{map}, either a string or the
#' column number. If \code{NULL} the first column is assumed to be the SNP name.
#' @param mapChromosome the column corresponding to the cromosome number in \code{map}, either a string or the
#' column number. If \code{NULL} the second column is assumed to be the chromosome number.
#' @param position the column corresponding to the position in \code{map}, either a string or the
#' column number. If \code{NULL} the third column is assumed to be the position.
#' @param kin string, specifying a csv file with the kinship matrix, with genotypes in both the rows
#' and the columns. Alternatively, an dataframe with a similar layout. Row names and column names should
#' be identical to column names in geno. If \code{NULL} a kinship matrix is computed from \code{geno}
#' using GRM.
#' @param pheno string, specifying a csv file with phenotypic data, with genotypes in both the rows
#' and the columns. Alternatively, an dataframe with a similar layout. Row names and column names should
#' be identical to column names in geno.
#' @param phenoGeno the column corresponding to the genotype in \code{pheno}, either a string or the
#' column number. If \code{NULL} the first column is assumed to be the genotype.
#'
#' @return an object of class gData with the following components:
#' \itemize{
#' \item{map a dataframe containing map data}
#' \item{markers a dataframe containing marker information}
#' \item{pheno a dataframe containing phenotypic data}
#' \item{kinship a kinship matrix}
#' \item{genotypes a charcter vector of genotypes}
#' \item{chromosomes a numeric vector of chromosome numbers}
#' \item{nChr the number of chromosomes}
#' \item{nGeno the number of genotypes}
#' \item{nSNP the number of SNPs}
#' \item{chrLengthsBp a numeric vector of chromosome lengths}
#' }

## TO DO: external
## TO DO: genes

gData <- function(geno,
  genoHeader = TRUE,
  genoRowNames = TRUE,
  map,
  mapSnpName = 1,
  mapChromosome = 2,
  mapPosition = 3,
  kin = NULL,
  pheno = NULL,
  phenoGeno = 1,
  createExternal = "none") {

  if (missing(geno) || (!is.data.frame(geno) && (!is.character(geno) && !file.exists(geno))))
    stop("geno should be either a valid file name or a dataframe")

  ## If geno is a character string, import corresponding file. Otherwise assign the data frame.
  if (is.character(geno)) {
    if (genoRowNames) rowNames = 1 else rowNames = NULL
    markers <- read.table(geno, sep = ",", header = genoHeader, row.names = rowNames)
  }
  else markers <- geno

  ## Check for row and column names. If not available give default names.
  if (length(rownames(markers)) == 0) {
    rownames(markers) <- paste0("g", formatC(1:nrow(markers),
      width = ceiling(log10(nrow(markers))), flag = "0"))}
  if (length(colnames(markers)) == 0) {
    colnames(markers) <- paste0("mk", formatC(1:ncol(markers),
      width = ceiling(log10(ncol(markers))), flag = "0"))}

  if (missing(map) || (!is.data.frame(map) && (!is.character(map) && !file.exists(map))))
    stop("map should be either a valid file name or a dataframe")
  if ((!is.numeric(mapSnpName) && !is.character(mapSnpName)) || length(mapSnpName) > 1 )
    stop("mapSnpName should be either an integer or a string")
  if ((!is.numeric(mapChromosome) && !is.character(mapChromosome)) || length(mapChromosome) > 1 )
    stop("mapChromosome should be either an integer or a string")
  if ((!is.numeric(mapPosition) && !is.character(mapPosition)) || length(mapPosition) > 1 )
    stop("mapPosition should be either an integer or a string")

  ## If map is a character string, import corresponding file. Otherwise assign the data frame.
  if (is.character(map)) {
    mapIn <- read.table(map, header = TRUE, stringsAsFactors = FALSE)
  } else {
    mapIn <- map}

  ## Select and rename the necessary columns from map.
  mapOut <- data.frame(snpName = mapIn[, mapSnpName],
    chromosome = mapIn[, mapChromosome],
    position = mapIn[, mapPosition],
    row.names = mapIn[, mapSnpName],
    stringsAsFactors = FALSE)

  ## Add cumulative position to map.
  chromosomes <- sort(unique(mapOut$chromosome))
  nChr <- length(chromosomes)
  chrLengths <- as.integer(table(mapOut$chromosome))

  if (nChr > 1) {
    chrPos <- c(0, cumsum(chrLengths)[1:(nChr - 1)])
    chrLengthsBp <- c(0, mapOut$position[cumsum(chrLengths)[1:(nChr - 1)]])
    cumPosition <- mapOut$position + rep(cumsum(chrLengthsBp), times = chrLengths)
  } else {
    cumPosition <- mapOut$position
    chrPos <- 0
    chrLengthsBp <- nrow(markers)
  }
  mapOut <- cbind(mapOut, cumPosition)
  mapOut <- mapOut[rownames(markers), ]

  ## If kin is a character string, import corresponding file. Otherwise assign the matrix.
  ## If kin is NULL compute the kinship matrix using GRM.
  if (!is.null(kin) && !is.matrix(kin) && (!is.character(kin) || !file.exists(kin)))
    stop("geno should be either NULL, a valid file name or a matrix")
  if (is.character(kin)) {
    kinship <- as.matrix(read.table(kin, sep = ",", header = TRUE))
    if (!identical(sort(rownames(kinship)), colnames(markers)) ||
      !identical(sort(colnames(kinship)), colnames(markers)))
      stop("row and column names of kin should be identical to row and column names of geno")
  } else {
    kinship <- GRM(t(markers))
  }

  ## If pheno is a character string, import corresponding file. Otherwise assign the data frame.
  if (!is.null(pheno) && !is.matrix(pheno) && (!is.character(pheno) || !file.exists(pheno)))
    stop("pheno should be either NULL, a valid file name or a data frame")
  if (is.character(pheno)) {
    phenoIn <- read.table(pheno, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    genoPos <- ifelse(is.character(phenoGeno), which(colnames(pheno) == phenoGeno), phenoGeno)
    colnames(pheno)[genoPos] <- "genotype"
  } else {
    phenoIn <- data.frame()
  }

  # if (!(createExternal %in% c('none','plink','scan_GLS'))) {createExternal <- 'none'}
  # if (createExternal == 'none')
  # {external <- list(bin.name="", kinship.name="", plink.name="")}
  # else if (createExternal == 'plink') {}
  # else if (createExternal == 'scan_GLS') {
  #   external <- list(bin.name = paste(description, ".bin", sep = ""),
  #     kinship.name = paste(description,".csv",sep=""), plink.name = '')
  #   MakeBin(kinship = kin, plant.names = names(markers), pheno = pheno,
  #     csv.file.name = geno,
  #     bin.file.name = external$bin.name)
  # }

  structure(list(map = mapOut,
    markers = markers,
    pheno = pheno,
    kinship = kinship,
    # external = external,
    genotypes = colnames(markers),
    chromosomes = chromosomes,
    nChr = nChr,
    nGeno = ncol(markers),
    nSNP = nrow(markers),
    chrLengthsBp = chrLengthsBp),
    # genes = data.frame()),
    class = "gData")
}



