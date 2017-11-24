#' S3 Class gData
#'
#' \code{createGData} creates an object of S3 class gData with genotypic and phenotypic data for usage in
#' further analysis. All input to the function is optional, however at least one input should be
#' provided. It is possible to provide an existing \code{gData} object as additional input in which case
#' data is added to this object. Existing data will be overwritten with a warning.\cr\cr
#' \code{is.gData} tests if an \code{R} object is of class \code{gData}.
#'
#' @param gData an optional gData object to be modified. If \code{NULL} a new gData object is created.
#' @param geno matrix or data.frame with genotypes in the rows and markers in the columns. A matrix from
#' the \code{matrix} in the base package may be provided as well as as matrix from the Matrix package. \cr
#' If no row names are used they are taken from \code{pheno} (if supplied and dimension matches).
#' If no column names are used the row names are taken from \code{map}
#' (if supplied and dimension matches).
#' @param map a data.frame with columns \code{chr} for chromosome and \code{pos} for position. Positions can
#' be in basepair or centimorgan. They should not be cumulative over the chromosomes.
#' Other columns are ignored. Marker names should be in the row names. These should match
#' the marker names in \code{geno} (if supplied).
#' @param kin a kinship matrix or list of kinship matrices with genotype in rows and colums.
#' These matrices can be from the \code{matrix} class as defined in the base package or
#' from the \code{dsyMatrix} class, the class of symmetric matrices in the Matrix package.\cr
#' The genotypes should be identical to the genotypes in \code{geno}.\cr
#' If a list of kinship matrices is provided these are supposed to be chromosome specific matrices.
#' In that case their names should match those of the names of the chromosomes in \code{map} or in
#' case no names are provided, the number of matrices should match the number of chromoses
#' in \code{map} in which case default names are provided.
#' @param pheno a data.frame or a list of data.frames with phenotypic data, with genotype in the
#' first column \code{genotype} and traits in the following columns. The trait columns should
#' be numerical columns only. A list of data.frames can be used for replications,
#' i.e. different environments.
#' @param covar a data.frame with extra covariates per genotype. Genotype should be in the rows.
#' @param x \code{R} object
#'
#' @return \code{createGData} returns an object of class \code{gData} with the following components:
#' \item{\code{map}}{a data.frame containing map data. Map is sorted by chromosome and position.}
#' \item{\code{markers}}{a sparse matrix from the Matrix package containing marker information in case
#' of numerical genotypic data, a standard matrix otherwise.}
#' \item{\code{pheno}}{a list of data.frames containing phenotypic data}
#' \item{\code{kinship}}{a kinship matrix of class \code{dsyMatrix} from the Matrix package.}
#' \item{\code{covar}}{a data.frame with extra covariates.}
#' \cr
#' \code{is.gData} returns \code{TRUE} or \code{FALSE} depending on whether its argument
#' is a \code{gData} object.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.gData}}
#'
#' @examples set.seed(1234)
#' ## Create genotypic data.
#' geno <- matrix(sample(x = c(0, 1, 2), size = 15, replace = TRUE), nrow = 3)
#' dimnames(geno) <- list(paste0("G", 1:3), paste0("M", 1:5))
#'
#' ## Construct map.
#' map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5, row.names = paste0("M", 1:5))
#'
#' ## Compute kinship matrix.
#' kin <- IBS(geno)
#'
#' ## Create phenotypic data.
#' pheno <- data.frame(paste0("G", 1:3), matrix(rnorm(n = 12, mean = 50, sd = 5), nrow = 3),
#'  stringsAsFactors = FALSE)
#' dimnames(pheno) = list(paste0("G", 1:3), c("genotype", paste0("T", 1:4)))
#'
#' ## Combine all data in gData object.
#' gData <- createGData(geno = geno, map = map, kin = kin, pheno = pheno)
#' summary(gData)
#'
#' ## Construct covariate.
#' covar <- data.frame(C1 = c("a", "a", "b"), row.names = paste0("G", 1:3))
#'
#' ## Compute alternative kinship matrix.
#' kin2 <- astle(geno)
#'
#' ## Add covariates to previously created gData object
#' gData2 <- createGData(gData = gData, kin = kin2, covar = covar)
#'
#' @name gData
NULL

#' @rdname gData
#' @export
createGData <- function(gData = NULL,
                        geno,
                        map,
                        kin,
                        pheno,
                        covar) {
  ## Check gData
  if (!is.null(gData) && !is.gData(gData)) {
    stop("Provided gData object should be of class gData.\n")
  }
  ## Check that at least one input argument, other than gData, is provided.
  if (missing(geno) && missing(map) && missing(kin) && missing(pheno) && missing(covar)) {
    stop("At least one of geno, map, kin, pheno and covar should be provided.\n")
  }
  ## Modify map
  if (!missing(map)) {
    if (!is.data.frame(map)) {
      stop("map should be a data.frame.\n")
    }
    if (!all(c("chr", "pos") %in% colnames(map))) {
      stop("chr and pos should be columns in map.\n")
    }
    if (!is.null(gData$map)) {
      warning("existing map will be overwritten.\n", call. = FALSE)
    }
    ## Extract columns and order
    map <- map[c("chr", "pos")]
    map <- map[order(map$chr, map$pos), ]
    if (all(rownames(map) == as.character(1:nrow(map)))) {
      ## If no marker name in input compute them from chromosome and position.
      ## Names are made unique if necessary by adding a suffix _1, _2, etc.
      replicates <- dplyr::count(map, map$chr, map$pos)$n
      suffix <- unlist(sapply(X = replicates, FUN = function(n) {
        if (n == 1) {
          return("")
        } else {
          return(paste0("_", 1:n))
        }
      }))
      rownames(map) <- paste0("chr", map$chr, "_", map$pos, suffix)
      warning("map contains no marker names. Names constructed from chromosome and position.\n",
              call. = FALSE)
    }
  } else {
    if (!is.null(gData$map)) {
      map <- gData$map
    } else {
      map <- NULL
    }
  }
  ## Modify pheno
  if (!missing(pheno)) {
    if (!is.data.frame(pheno) &&
        !(is.list(pheno) &&
          all(sapply(X = pheno, FUN = is.data.frame)))) {
      stop("pheno should be a data.frame or a list data.frames.\n")
    }
    if (is.data.frame(pheno)) {
      ## If not a list already put data.frame/matrix in a list.
      pheno <- setNames(list(pheno), deparse(substitute(pheno)))
    } else {
      if (is.null(names(pheno))) {
        ## Add default names.
        names(pheno) <- sapply(X = 1:length(pheno), FUN = function(x) {
          paste0("Environment", x)
          })
        message("pheno contains no environment names. Default names added.\n")
      } else {
        if (!isTRUE(all(sapply(X = names(pheno), FUN = nchar) > 0))) {
          ## Add default names for unnamed environments.
          names(pheno) <- sapply(X = 1:length(pheno), FUN = function(x) {
            if (!isTRUE(nchar(names(pheno)[x]) > 0)) {
              paste0("Environment", x)
            } else {
              names(pheno)[x]
            }
          })
          message("Some data.frames in pheno contain no environment names.
                  Default names added.\n")
        }
      }
    }
    if (!all(sapply(X = pheno, FUN = function(x) {
      colnames(x)[1] == "genotype"
      }))) {
      stop("First column in pheno should be genotype.\n")
    }
    ## Convert genotype to character.
    for (i in 1:length(pheno)) {
      pheno[[i]]$genotype <- as.character(pheno[[i]]$genotype)
    }
    ## Check that all non-genotype columns are numerical.
    if (!all(sapply(X = pheno, FUN = function(x) {
      all(sapply(X = x[-1], FUN = is.numeric))
    }))) {
      stop("all trait columns in pheno should be numerical.\n")
    }
    if (!is.null(gData$pheno)) {
      warning("existing pheno will be overwritten.\n", call. = FALSE)
    }
  } else {
    if (!is.null(gData$pheno)) {
      pheno <- gData$pheno
    } else {
      pheno <- NULL
    }
  }
  ## Modify geno
  if (!missing(geno)) {
    if (!is.data.frame(geno) && !inherits(geno, "Matrix") && !is.matrix(geno)) {
      stop("geno should be a matrix or a data.frame.\n")
    }
    if (is.data.frame(geno) || is.matrix(geno)) {
      if (is.numeric(unlist(geno))) {
      ## Convert geno to Matrix of class Matrix.
      markers <- as(geno, "Matrix")
      } else {
        markers <- as.matrix(geno)
      }
    } else {
      markers <- geno
    }
    ## Check for row names in markers. If not available take them from pheno or use default names.
    if (all(rownames(markers) == as.character(1:nrow(markers)))) {
      if (missing(pheno) || is.null(pheno)) {
        ## Default names are constructed as g001, g002, etc. with the number of 0 dependent on the
        ## number of rows.
        rownames(markers) <- paste0("g",
                                    formatC(1:nrow(markers),
                                            width = ceiling(log10(nrow(markers))),
                                            flag = "0")
                                    )
        warning("geno contains no genotype names. Default names used.\n", call. = FALSE)
      } else {
        if (nrow(pheno[[1]]) == nrow(markers)) {
          rownames(markers) <- rownames(pheno[[1]])
          warning("geno contains no genotype names. Names taken from pheno.\n", call. = FALSE)
        } else {
          stop("geno contains no genotype names. Dimensions between geno and pheno differ.\n")
        }
      }
    } else {
      ## Sort rownames alphabetically.
      markers <- markers[order(rownames(markers)), , drop = FALSE]
    }
    if (is.null(colnames(markers))) {
      ## Check for column names in markers. If not available take them from map.
      if (is.null(map)) {
        stop("geno contains no marker names. Map not available.\n")
      }
      if (nrow(map) != ncol(markers)) {
        stop("geno contains no marker names. Dimensions between geno and map differ.\n")
      }
      colnames(markers) <- rownames(map)
      warning("geno contains no marker names. Names taken from map.\n", call. = FALSE)
    } else if (!is.null(map)) {
      if (any(!colnames(markers) %in% rownames(map[rownames(map) %in% colnames(markers), ]))) {
        warning("not all markers in geno are in map. Extra markers are removed.\n", call. = FALSE)
      }
      markers <- markers[, colnames(markers) %in% rownames(map)]
    }
    if (!is.null(gData$markers)) {
      warning("existing geno will be overwritten.\n", call. = FALSE)
    }
  } else {
    if (!is.null(gData$markers)) {
    markers <- gData$markers
    } else {
      markers <- NULL
    }
  }
  ## Modify kin
  if (!missing(kin)) {
    if (!is.null(kin) &&
        !inherits(kin, "Matrix") &&
        !is.matrix(kin) &&
        !(is.list(kin) &&
          all(sapply(X = kin, FUN = function(x) {
          inherits(x, "Matrix") || is.matrix(x)
        })))) {
      stop("kin should be a matrix or a list of matrices.\n")
    }
    if (!is.null(map) && is.list(kin) && length(kin) != dplyr::n_distinct(map$chr)) {
      stop("kin should be the same length as the number of chromosomes in map.\n")
    }
    if (is.null(rownames(kin)) || is.null(colnames(kin))) {
      stop("row and column names in kin cannot be NULL.\n")
    }
    if ((!is.null(markers) && is.list(kin) &&
        any(sapply(X = kin, FUN = function(x) {
          !all(rownames(x) %in% rownames(markers)) || !all(colnames(x) %in% rownames(markers))
        }))) ||
        (!is.null(markers) && !is.list(kin) &&
         (!all(rownames(kin) %in% rownames(markers)) ||
          !all(colnames(kin) %in% rownames(markers))))) {
      stop("row and column names of kin should be in row and column names of geno.\n")
    }
    ## Order as in geno.
    if (!is.null(names(kin)) && names(kin) != unique(map$chr)) {
      stop("names of kin should correspond to names of chromosomes in map.\n")
    }
    if (is.list(kin)) {
      kin <- lapply(X = kin,
                    FUN = function(x) {
                      as(x[order(match(rownames(x), rownames(markers))),
                           order(match(colnames(x), rownames(markers)))],
                         "dsyMatrix")
                    })
    } else {
      if (is.matrix(kin)) {
        kin <- as(kin, "dsyMatrix")
      }
      kin <- kin[order(match(rownames(kin), rownames(markers))),
                 order(match(colnames(kin), rownames(markers)))]
    }
    ## Add default names.
    if (is.list(kin) && is.null(names(kin))) {
      warning("kin contains no names. Default names added.\n")
      names(kin) <- unique(map$chr)
    }
    if (!is.null(gData$kinship)) {
      warning("existing kinship will be overwritten.\n", call. = FALSE)
    }
  } else {
    if (!is.null(gData$kinship)) {
      kin <- gData$kinship
    } else {
      kin <- NULL
    }
  }
  ## Modify covar
  if (!missing(covar)) {
    if (!is.null(covar) && !is.data.frame(covar)) {
      stop("covar should be a data.frame.\n")
    }
    if (!all(sapply(X = covar, FUN = function(x) {
      is.numeric(x) || is.character(x) || is.factor(x)}))) {
      stop("all columns in covar should be numeric, character or factor columns.\n")
    }
    ## Convert character columns to factors.
    covar[sapply(covar, is.character)] <- lapply(X = covar[sapply(X = covar,
                                                                  FUN = is.character)],
                                                 FUN = as.factor)
    if (!is.null(gData$covar)) {
      warning("existing covar will be overwritten.\n", call. = FALSE)
    }
  } else {
    if (!is.null(gData$covar)) {
      covar <- gData$covar
    } else {
      covar <- NULL
    }
  }
  ## Create gData object.
  gData <- structure(list(map = map,
                          markers = markers,
                          pheno = pheno,
                          kinship = kin,
                          covar = covar),
                     class = "gData")
  return(gData)
}

#' @rdname gData
#' @export
is.gData <- function(x) {
  inherits(x, "gData")
}

