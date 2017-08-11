#' Create a gData object
#'
#' \code{createGData} creates an object of S3 class with genotypic and phenotypic data for usage in
#' further analysis. All input to the function is optional, however at least one input should be
#' provided.
#'
#' @param geno matrix or data.frame with genotypes in the rows and markers in the columns.
#' If no row names are used they are taken from \code{pheno} (if supplied and dimension matches).
#' If no column names are used the row names are taken from \code{map}
#' (if supplied and dimension matches).
#' @param map a data.frame with columns chr for chromosome and pos for position. Positions can
#' be in basepair or centimorgan. They should not be cumulative over the chromosomes.
#' Other columns are ignored. Marker names should be in the row names. These should match
#' the marker names in \code{geno} (if supplied).
#' @param kin a kinship matrix with genotype in rows and colums. These should be identical to the genotypes
#' in \code{geno}
#' @param pheno a data.frame or a list of data.frames with phenotypic data, with genotype in the
#' first column and traits in the following columns. A list of data.frames can be used for replications,
#' i.e. different environments.
#' @param covar a data.frame with extra covariates per genotype. Genotype should be in the rows.
#'
#' @return an object of class gData with the following components:
#' \itemize{
#' \item{map a data.frame containing map data. Map is sorted by chromosome and position.}
#' \item{markers a matrix containing marker information}
#' \item{pheno a list of matrices containing phenotypic data}
#' \item{kinship a kinship matrix}
#' \item{covar a data.frame with extra covariates}
#' }

createGData <- function(geno,
  map,
  kin,
  pheno,
  covar) {
  ## Check that at least one input is provided.
  if(missing(geno) && missing(map) && missing(kin) && missing(pheno) && missing(covar))
    stop("At least one of geno, map, kin, pheno and covar should be provided.\n")
  ## Modify map
  if (!missing(map)) {
    if (!is.data.frame(map)) stop("map should be a data.frame.\n")
    if (!all(c("chr", "pos") %in% colnames(map))) stop("chr and pos should be columns in map.\n")
    ## Extract columns and order
    map <- map[c("chr", "pos")]
    map <- map[order(map$chr, map$pos), ]
    if (all(rownames(map) == as.character(1:nrow(map)))) {
      ## If no marker name in input compute them from chromosome and position.
      ## Names are made unique if necessary by adding a suffix _1, _2, etc.
      replicates <- dplyr::count(map, chr, pos)$n
      suffix <- unlist(sapply(replicates, FUN = function(n) {
        if (n == 1) return("") else return(paste0("_", seq(1:n)))}))
      rownames(map) <- paste0("chr", map$chr, "_", map$pos, suffix)
      warning("map contains no marker names. Names constructed from chromosome and position.\n",
        call. = FALSE)
    }
  } else map <- NULL
  ## Modify pheno
  if (!missing(pheno)) {
    if (!is.data.frame(pheno) &&
        !(is.list(pheno) && all(sapply(pheno, FUN = is.data.frame))))
      stop("pheno should be a data.frame or a list data.frames.\n")
    if (is.data.frame(pheno)) {
      ## If not a list already put data.frame/matrix in a list.
      pheno <- list(pheno)
    }
    if (!all(sapply(pheno, FUN = function(x) {colnames(x)[1] == "genotype"})))
      stop("First column in pheno should be genotype.\n")
    ## Convert all list items to matrices.
  } else pheno <- NULL
  ## Modify geno
  if (!missing(geno)) {
    if (!is.data.frame(geno) && !is.matrix(geno)) stop("geno should be a matrix or a data.frame.\n")
    if (is.data.frame(geno)) markers <- as.matrix(geno) else markers <- geno
    ## Check for row names in markers. If not available take them from pheno or use default names.
    if (all(rownames(markers) == as.character(1:nrow(markers)))) {
      if (missing(pheno)) {
        ## Default names are constructed as g001, g002, etc. with the number of 0 dependent on the
        ## number of rows.
        rownames(markers) <- paste0("g", formatC(1:nrow(markers),
          width = ceiling(log10(nrow(markers))), flag = "0"))
        warning("geno contains no genotype names. Default names used.\n", call. = FALSE)
      } else if (nrow(pheno[[1]]) == nrow(markers)) {
        rownames(markers) <- rownames(pheno[[1]])
        warning("geno contains no genotype names. Names taken from pheno.\n", call. = FALSE)
      } else {
        stop("geno contains no genotype names. Dimensions between geno and pheno differ.\n")
      }
    } else {
      ## Sort rownames alphabetically.
      markers <- markers[order(rownames(markers)), ]
    }
    if (is.null(colnames(markers))) {
      ## Check for column names in markers. If not available take them from map.
      if (missing(map)) stop("geno contains no marker names. Map not available.\n")
      if (nrow(map) != ncol(markers))
        stop("geno contains no marker names. Dimensions between geno and map differ.\n")
      colnames(markers) <- rownames(map)
      warning("geno contains no marker names. Names taken from map.\n", call. = FALSE)
    }
  } else markers <- NULL
  ## Modify kin
  if (!missing(kin)) {
    if (!is.null(kin) && !is.matrix(kin))
      stop("kin should be a matrix.\n")
    if (!missing(geno) &&
        (!all(rownames(kin) %in% rownames(markers)) ||
            !all(colnames(kin) %in% rownames(markers))))
      stop("row and column names of kin should be in row and column names of geno.\n")
    ## Order as in geno.
    kin <- kin[order(match(rownames(kin), rownames(markers))), order(match(colnames(kin), rownames(markers)))]
  } else kin <- NULL
  ## Modify covar
  if (!missing(covar)) {
    if (!is.null(covar) && !is.data.frame(covar))
      stop("covar should be a data.frame.\n")
  } else covar <- NULL
  ## Create gData object.
  gData <- structure(list(map = map,
    markers = markers,
    pheno = pheno,
    kinship = kin,
    covar = covar),
    class = "gData")
  return(gData)
}

is.gData <- function(x) {
  inherits(x, "gData")
}



