#' Create a gData object
#'
#' Create an object of S3 class gData based on external files or dataframes containing genotypic and phenotypic
#' information.
#'
#' Using the argument \code{pheno} can only be used for adding phenotypic data in a dataframe. For adding
#' phenotypic data from a file use \code{\link{addPhenoData}}.
#'
#' @param geno string, specifying a csv file with the genotypic data, genotypes in the rows and
#' markers in the columns. Alternatively, a dataframe with a similar layout.
#' @param genoHeader does \code{geno} contain a header row. If \code{FALSE} markers are given a default
#' name derived from the \code{map}.
#' @param genoRowNames does \code{geno} contain rownames in its first column. If \code{FALSE} genotypes
#' are given a default name.
#' @param map string, specifying a csv file with at least columns for chromosome and position. Column for
#' SNP Name can be provided as well. The positions should be in base-pair or centimorgan. They should not
#' be cumulative over the chromosomes. Other columns are ignored. Alternatively, a dataframe with a
#' similar layout.
#' @param mapSnpName the column corresponding to the SNP name in \code{map}, either a string or the
#' column number. If \code{NULL} the SNP name is constructed from the chromosome and position.
#' @param mapChromosome the column corresponding to the cromosome number in \code{map}, either a string or the
#' column number. Default the second column is assumed to be the chromosome number.
#' @param mapPosition the column corresponding to the position in \code{map}, either a string or the
#' column number. Default the third column is assumed to be the position.
#' @param kin string, specifying a csv file with the kinship matrix, with genotypes in both the rows
#' and the columns. Alternatively, a matrix with a similar layout. Row names and column names should
#' be identical to column names in geno. If \code{NULL} a kinship matrix is computed from \code{geno}
#' using GRM.
#' @param pheno a dataframe with phenotypic data, with genotypes in the rows and traits in the columns.
#' Genotypes in pheno should be in column names in geno. See details.
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
#' \item{chrLengths a numeric vector of chromosome lengths}
#' }

gData <- function(geno,
  genoHeader = TRUE,
  genoRowNames = TRUE,
  map,
  mapSnpName = NULL,
  mapChromosome = 2,
  mapPosition = 3,
  kin = NULL,
  pheno = NULL,
  phenoGeno = 1) {

  if (missing(geno) || (!is.data.frame(geno) && (!is.character(geno) && !file.exists(geno))))
    stop("geno should be either a valid file name or a dataframe\n")

  ## If geno is a character string, import corresponding file. Otherwise assign the data frame.
  if (is.character(geno)) {
    if (genoRowNames) rowNames = 1 else rowNames = NULL
    markers <- read.table(file = geno, header = genoHeader, sep = ",", row.names = rowNames)
  }
  else markers <- geno

  if (missing(map) || (!is.data.frame(map) && (!is.character(map) && !file.exists(map))))
    stop("map should be either a valid file name or a dataframe\n")
  if ((!is.null(mapSnpName) && !is.numeric(mapSnpName) && !is.character(mapSnpName)) || length(mapSnpName) > 1)
    stop("mapSnpName should be either NULL, an integer or a string\n")
  if ((!is.numeric(mapChromosome) && !is.character(mapChromosome)) || length(mapChromosome) > 1)
    stop("mapChromosome should be either an integer or a string\n")
  if ((!is.numeric(mapPosition) && !is.character(mapPosition)) || length(mapPosition) > 1)
    stop("mapPosition should be either an integer or a string\n")

  ## If map is a character string, import corresponding file. Otherwise assign the data frame.
  if (is.character(map)) {
    mapIn <- read.table(file = map, header = TRUE, stringsAsFactors = FALSE)
  } else {
    mapIn <- map}

  ## Check that mapChromosome and mapPosition are columns in map
  if ((is.character(mapChromosome) && !mapChromosome %in% colnames(mapIn)) ||
      (is.numeric(mapChromosome) && mapChromosome > ncol(mapIn)) ||
      any(mapIn[, mapChromosome] != round(mapIn[, mapChromosome])))
    stop("mapChromosome should an integer column in map")
  if ((is.character(mapPosition) && !mapPosition %in% colnames(mapIn)) ||
      (is.numeric(mapPosition) && mapPosition > ncol(mapIn)))
    stop("mapPosition should a column in map")

  ## Select and rename the necessary columns from map.
  mapOut <- data.frame(
    chromosome = mapIn[, mapChromosome],
    position = mapIn[, mapPosition],
    stringsAsFactors = FALSE)

  mapOut <- dplyr::arrange(mapOut, chromosome, position)
  ## If no snpName in input, extract it from markers and otherwise compute it from chromosome and position.
  if (is.null(mapSnpName)) {
    if (genoHeader) {
      if (nrow(mapOut) != ncol(markers))
        stop("map contains no SNP names. Number of SNPs in map should be identical to number of
          SNPs in geno.\n")
      mapOut <- tibble::add_column(mapOut, snpName = colnames(markers), .before = 1)
      warning("map contains no SNP names. Names taken from markers.\n", call. = FALSE)
    } else {
      replicates <- dplyr::count(mapOut, chromosome, position)$n
      suffix <- unlist(sapply(replicates, FUN = function(n) {
        if (n == 1) return("") else return(paste0("_", seq(1:n)))}))
      mapOut <- tibble::add_column(mapOut,
        snpName = paste0("chr", mapOut$chromosome, "_", mapOut$position, suffix),
        .before = 1)
      warning("map contains no SNP names. Names composed from chromosome and position.\n", call. = FALSE)
    }
  } else {
    ## Check that mapSnpName is column in map
    if ((is.character(mapSnpName) && !mapSnpName %in% colnames(mapIn)) ||
        (is.numeric(mapSnpName) && mapSnpName > ncol(mapIn)))
      stop("mapSnpName should a column in map.\n")
    mapOut <- tibble::add_column(mapOut, snpName = mapIn[, mapSnpName], .before = 1)
  }
  rownames(mapOut) <- mapOut$snpName

  ## Check for row and column names in markers. If not available give default names.
  ## Use genoRowNames for check on row names since default row names are always used.
  if (!genoRowNames) {
    rownames(markers) <- paste0("g", formatC(1:nrow(markers),
      width = ceiling(log10(nrow(markers))), flag = "0"))
    warning("geno contains no genotype names. Default names used.\n", call. = FALSE)
  }
  ## Use genoHeader for check on row names since default column names are always used.
  if (!genoHeader) {
    if (nrow(mapOut) != ncol(markers))
      stop("geno contains no SNP names. Number of SNPs in geno should be identical to number of
        SNPs in map.\n")
    colnames(markers) <- rownames(mapOut)
    warning("geno contains no SNP names. Names taken from map.\n", call. = FALSE)
  }
  ## Map may only contain SNPs which are also in markers.
  mapOut <- mapOut[which(mapOut$snpName %in% colnames(markers)), ]

  ## Add cumulative position to map.
  chromosomes <- sort(unique(mapOut$chromosome))
  nChr <- length(chromosomes)
  if (nChr > 1) {
    chrLengths <- c(0, as.numeric(aggregate(mapOut$position, by = list(mapOut$chromosome),
      FUN = max)$x))[1:nChr]
    cumPosition <- mapOut$position + rep(cumsum(chrLengths),
      times = as.integer(table(mapOut$chromosome)))
  } else {
    cumPosition <- mapOut$position
    chrLengths <- ncol(markers)
  }
  mapOut <- tibble::add_column(mapOut, cumPosition)

  ## If kin is a character string, import corresponding file. Otherwise assign the matrix.
  ## If kin is NULL compute the kinship matrix using GRM.
  if (!is.null(kin) && !is.matrix(kin) && (!is.character(kin) || !file.exists(kin)))
    stop("geno should be either NULL, a valid file name or a matrix\n")
  if (is.character(kin)) {
    kinship <- as.matrix(read.table(file = kin, header = TRUE, sep = ","))
    if (!identical(sort(rownames(kinship)), rownames(markers)) ||
        !identical(sort(colnames(kinship)), rownames(markers)))
      stop("row and column names of kin should be identical to row and column names of geno.\n")
  } else {
    kinship <- GRM(markers)
  }

  ## Assign the phenotypic data frame to pheno.
  if (!is.null(pheno) && !is.data.frame(pheno))
    stop("pheno should be either NULL or a data frame.\n")
  if (is.data.frame(pheno)) {
    genoPos <- ifelse(is.character(phenoGeno), which(colnames(pheno) == phenoGeno), phenoGeno)
    colnames(pheno)[genoPos] <- "genotype"
    if (!all(pheno$genotype) %in% rownames(markers))
      stop("all genotypes in pheno should be in geno.\n")
  } else {
    pheno <- data.frame()
  }

  structure(list(map = mapOut,
    markers = markers,
    pheno = pheno,
    kinship = kinship,
    genotypes = rownames(markers),
    chromosomes = chromosomes,
    nChr = nChr,
    nGeno = ncol(markers),
    nSNP = nrow(markers),
    chrLengths = chrLengths),
    class = "gData")
}

is.gData <- function(x) {
  inherits(x, "gData")
}



