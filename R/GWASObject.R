#'  OBJECTIVE : given several internal and/or external objects, create a GWAS object
#
#' @param markers a string, specifying a csv-file name with the markerdata (row names: marker names;
#' column names: genotype names). Alternatively, an dataframe with similar layout.
#' @param map a dataframe with two columns, named chromosome and position.
#' The positions should be in base-pair or morgan. They should not be cumulative over the chromosomes.
#' If NULL it is assumed the map information is contained in the \code{marker} data in columns named
#' chromosome and position.

#' @param  description: name of the data-set
#' @param  pheno : an R-data-frame with the phenotypic data. Importing pheno-data from a text/csv file can be done
#'         using the AddPhenoData function
#' @param  accession.name.type : the type of accession names used in marker.object:
#'    1= "stock.new", 2= "stock.old", 3= "ecotype.standard", 4="ecotype.yi", 5="array.id"
#' @param  allele.name.type   : the allele-encoding used in marker.object:
#'    1="columbia as one",2="minor allele as one",3="other"
#' @param  marker.name.type   : the type of marker names used in marker.object:
#'    1="chromosome, position and gene",2="numbers",3="other"
#' @param maf : minor allele frequency to which the markers in marker.object have been restricted
#'          (in particular : before, when using the ReadMarkerAndGeneData function)
#
# * csvName
# * binName
# * RimageName : name of the R-image the GWASobject will be saved in
# * returnObject : if true, the GWASobject is returned, otherwise it is only written to the RData file
#
# OUTPUT :
# * a GWAS object (if returnObject=TRUE), where many entries and filenames will start, with "description "
# *
#
################################


GWASObject <- function(markers = NULL, map = NULL, kinship = NULL, pheno = NULL,
  kinshipType = "grm") {

  if (!is.character(markers) && !is.data.frame(markers))
    stop("markers should be either a string or a dataframe")

  if (is.character(markers)) {
    markers <- read.table(markers)
    ## Check for row and column names. If not available give default names.
  } else {
    if (length(rownames(markers)) == 0)
    {rownames(markers) <- paste0("ind", 1:nrow(markers))}
    if (length(colnames(markers)) == 0)
    {names(markers) <- paste0("m", 1:ncol(markers))}
  }

  if (all(c("chromosome", "position") %in% colnames(markers))) {
    map <- markers[c("chromosome", "position")]
    markers <- markers[, -which(colnames(markers) %in% c("chromosome", "position"))]
  } else {
    if (!is.null(map) && is.data.frame(map) && all(c("chromosome","position") %in% names(map))) {
      map <- map[c("chromosome", "position")]
    } else {
      map <- data.frame(chromosome = rep(1, nrow(markers)), position = 1:nrow(markers))
    }
  }

  n <- ncol(markers)

  if (!is.null(kinship) && ncol(kinship) == n) {
    K <- kinship
  } else {
    if (!kinshipType %in% c("grm", "ibs", "id")) kinshipType <- "grm"
    if (kinshipType == "grm") {
      K <- GRM(t(markers))
    } else if (kinshipType == "ibs") {
      K <- IBS(markers)
    }
    if (kinshipType == "id")  {
      K <- diag(n)
    }
  }

  structure(
    list(markers = markers,
      map = map,
      K = K),
    class = "GWASObject"
  )
}



# MakeGwasObject <- function(marker.object,description="test",kinship=matrix(0),map=data.frame(),
#   pheno=data.frame(),gene.info=character(),kinship.type="ibs",
#   accession.name.type=1,allele.name.type=1,marker.name.type=2,maf=0,
#   gene.dataframe=data.frame(),generate.csv=FALSE,generate.bin=FALSE,
#   csvName=paste(description,".csv",sep=""),binName=paste(description,".bin",sep=""),
#   RimageName=paste(description,".RData",sep=""),returnObject=TRUE) {
#
#
#
#
#
#
#
#   ###
#
#   nchr            <- length(unique(map$chromosome))
#   chromosomes     <- unique(map$chromosome)
#   chr.lengths     <- as.numeric(table(map$chromosome))
#   if (nchr > 1) {
#     chr.pos         <- c(0,cumsum(chr.lengths)[1:(nchr-1)])
#     chr.lengths.bp  <- c(0,map$position[cumsum(chr.lengths)[1:(nchr-1)]])
#     cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
#   } else {
#     cumpositions    <- map$position
#     chr.pos         <- 0
#     chr.lengths.bp  <- nrow(marker.object)
#   }
#   plant.names <- names(marker.object)
#   N           <- nrow(marker.object)
#
#   ### principal components, as in Patterson et al 2006
#
#   gc()
#   PCAmatrix  <- ComputePcaMatrix(marker.object=marker.object,maf=maf)
#   # requiring too much memory :
#   # PCAmatrix  <- ComputePcaMatrix2(marker.object=marker.object,maf=maf)
#   PCAs       <- as.data.frame(ComputePcas(PCAmatrix)$pcas)
#   row.names(PCAs) <- plant.names
#   names(PCAs)<- paste("pca",as.character(1:ncol(PCAs)),sep="")
#   #PCAs <- data.frame()
#   ### kinship matrix
#
#   K.name          <- paste(description,"_","kinship.csv",sep="")
#   write.table(A, file=K.name,quote = FALSE, sep = ",", row.names = FALSE,col.names = plant.names)
#   AINV            <- MakeKinshipAsreml(A,genotype.names=plant.names)
#
#   ### create mapframe
#
#   if (length(gene.info)>0) {
#     map.frame <- cbind(map=map,data.frame(cum.position=cumpositions),gene1=gene.info$gene1,gene2=gene.info$gene2)
#   } else {
#     map.frame <- cbind(map=map,data.frame(cum.position=cumpositions))
#   }
#
#   ###
#
#   csv.name  <- csvName
#   bin.name  <- binName
#
#   if (generate.csv) {
#     csv.name=csvName
#     MakeCsv(markers=marker.object,plant.names=plant.names,file.name=csv.name)
#   } else {
#     csv.name=""
#   }
#   if (generate.bin) {
#     bin.name=binName
#     #MakeBin(kinship.file=K.name,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
#     #MakeBin(kinship=K.name,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
#     MakeBin(kinship=A,plant.names=plant.names,csv.file.name=csv.name,bin.file.name=bin.name) # old argument : markers=marker.object,
#   } else {
#     bin.name=""
#   }
#
#   # Define various factors
#   accession.levels   <- c("stock.new","stock.old","ecotype.standard","ecotype.yi","array.id")
#   allele.levels      <- c("columbia as one","minor allele as one","other")
#   marker.name.levels <- c("chromosome, position and gene","numbers","other")
#
#   GWAS.obj  <- list(description=description,markers=marker.object,pheno=pheno,map=map.frame,
#     kinship=A,kinship.asreml=AINV,genes=gene.dataframe,plant.names=plant.names,N=N,
#     n=n,nchr=nchr,chromosomes=chromosomes,chr.lengths.bp=chr.lengths.bp,
#     external=list(kinship.name=K.name,csv.name=csv.name,bin.name=bin.name),pca=PCAs,
#     markerInfo=list(accession.name.type=factor(x = accession.levels[accession.name.type], levels=accession.levels),
#       allele.name.type=factor(x = allele.levels[allele.name.type], levels=allele.levels),
#       marker.name.type=factor(x = marker.name.levels[marker.name.type], levels=marker.name.levels),maf=maf),
#     real.effects=list(locations=data.frame(NULL),sizes=data.frame(NULL))
#   )
#
#   names(GWAS.obj$map)[1:2] <- c("chromosome","position")
#
#   if (RimageName!="") {save(GWAS.obj,file=RimageName)}
#
#   if (returnObject) {return(GWAS.obj)}
#
# }



















is.GWASObject <- function(x) {
  inherits(x, "GWASObject")
}
