#' Add phenotypic data to a gData object
#'
#' Add phenotypic data from a csv file to a \code{\link{gData}} Object. The csv file should contain the genotype as
#' its first column. Other columns are imported as numerical except for those in \code{factorCols} which are
#' imported as factors.
#'
#' @param gData an object of class \code{\link{gData}} to which the phenotypic data should be added.
#' @param pheno string, specifying a csv file with genotype in the first column.
#' @param factorCols a vector of columns that should be imported as factors from the csv file. Either a
#' vector of column names or of column numbers.
#'
#' @return the input object of class \code{\link{gData}} with the imported phenotypic data in the
#' component pheno.

addPhenoData  <- function(gData,
  pheno,
  factorCols = NULL) {

  if (missing(gData) || !is.gData(gData))
    stop("gData should be a gData object\n")

  if (missing(pheno) || !is.character(pheno) || !file.exists(pheno))
    stop("pheno should be a valid file name\n")

  if (!is.null(factorCols) && !is.character(factorCols) && !is.numeric(factorCols))
    stop("factorCols should be either NULL, a vector of strings or a vector of integers\n")

  ## Read only the first line of data to get the column names.
  phenoIn <- read.table(file = pheno, header = TRUE, sep = ",", nrows = 1)
  ## Set the column classes.
  colClasses <- c("character", rep("numeric", ncol(phenoIn) - 1))

  if (!is.null(factorCols)) {
    if (is.character(factorCols) && !all(factorCols %in% colnames(phenoIn))) {
      stop("all factorCols should be column names in pheno")
    } else if (is.numeric(factorCols) && any(factorCols > ncol(phenoIn))) {
      stop("all factorCols should be columns in pheno")
    }
    factorPos <- ifelse(is.character(factorCols), which(colnames(phenoIn) %in% factorCols), factorCols)
    colClasses[factorPos] <- "factor"
  }
  ## Read the full table
  phenoIn <- read.table(file = pheno, sep=",", na.strings = "", header = TRUE, colClasses = colClasses)
  names(phenoIn)[1] <- "genotype"

  ## Remove the accessions which are not in genotypes
  excludeGenotypes <- unique(phenoIn$genotype[which(!phenoIn$genotype %in% gData$genotypes)])
  phenoOut <- phenoIn[which(phenoIn$genotype %in% gData$genotypes), ]

  ## Order by order in genotypes
  dplyr::right_join(phenoOut, data.frame(genotype = gData$genotypes), by = "genotype")

  gData$pheno <- phenoOut

  return(gData)
}
