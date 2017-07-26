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

## TO DO: Complete the genotype names?
## TO DO: addVarMeans
## TO DO: addNormalTransform
## TO DO: missings in factor columns

addPhenoData  <- function(gData,
  pheno,
  addVarMeans = FALSE,
  meanCols = 0,
  addNormalTransform = FALSE,
  factorCols = NULL) {

  if (missing(gData) || !is.gData(gData))
    stop("gData should be a gData object\n")

  if (missing(pheno) || !is.character(pheno) || !file.exists(pheno))
    stop("pheno should be a valid file name\n")

  if (!is.null(factorCols) && !is.character(factorCols) && !is.numeric(factorCols))
    stop("factorCols should be either NULL, a vector of strings or a vector of integers\n")

  if (!addVarMeans) {meanCols <- 0}
  if (sum(meanCols)==0) {addVarMeans<-FALSE}

  ## Read only the first line of data to get the column names.
  phenoIn <- read.table(file = pheno, header = TRUE, sep = ",", nrows = 1)
  ## Set the column classes.
  colClasses <- c("character", rep("numeric", ncol(phenoIn) - 1))

  if (!is.null(factorCols)) {
    if (is.character(factorCols) && !all(factorCols %in% colnames(phenoIn))) {
      stop("all factorCols should be column names in pheno")
    } else if (is.numeric(factorCols) && !all(factorCols <= ncol(phenoIn))) {
      stop("all factorCols should be columns in pheno")
    }
    factorPos <- ifelse(is.character(factorCols), which(colnames(phenoIn) %in% factorCols), factorCols)
    colClasses[factorPos] <- "factor"
  }
  ## Read the full table
  phenoIn <- read.table(file = pheno, sep=",", na.strings = "", header = TRUE, colClasses = colClasses)

  ## Complete the genotype names
  for (i in 2:nrow(phenoIn)) {
    if (is.na(phenoIn[i, 1])) {phenoIn[i, 1] <- phenoIn[i - 1, 1]}}
  names(phenoIn)[1] <- "genotype"

  ## Remove the accessions which are not in genotypes
  phenoOut <- phenoIn[which(phenoIn$genotype %in% gData$genotypes), ]

  ## Add extra accessions, i.e. the ones that do not occur in phenoOut, but DO occur in genotypes
  extraAcc <- setdiff(gData$genotypes, phenoOut$genotype)
  if (length(extraAcc) > 0) {
    phenoOut <- dplyr::add_row(phenoOut, genotype = extraAcc)
  }

  ## Sort to match sorting done automatically by table command
  phenoOut <- dplyr::arrange(phenoOut, by = genotype)

  ## Count the number of replicates in the remaining genotypes
  ## (regardless of missing values in the phenotype)
  replicates <- as.numeric(table(phenoOut$genotype))

  ## Add _1, _2, etc to the row names to make them unique over replicates
  row.names(phenoOut) <- paste0(phenoOut$genotype, "_" ,
    unlist(sapply(replicates, FUN = function(n) {seq(1:n)})))

  ## Reorder by order in genotypes
  dplyr::right_join(phenoOut, data.frame(genotype = replicates), by = "genotype")



  if (addVarMeans) {
    NC    <- ncol(phenoOut)
    phenoOut <- AddMeans(input.frame=phenoOut,col.select=meanCols)   # the averages of the columns given in col.select are now added as extra columns
    if (addNormalTransform) {
      for (i in c(meanCols,NC+1:length(meanCols))) {
        phenoOut <- cbind(phenoOut,qqnorm(phenoOut[,i],plot.it=F)$x)
        names(phenoOut)[ncol(phenoOut)] <- paste(names(phenoOut)[i],"_transformed",sep="")
      }
    }
  }

  if (!is.null(factorCols)) {
    for (i in factorCols) {
      phenoOut[is.na(phenoOut[, i]), i] <- levels(phenoOut[, i])[1]
    }
  }

  gData$pheno <- phenoOut

  return(gData)
}
