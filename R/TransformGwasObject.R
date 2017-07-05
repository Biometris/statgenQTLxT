TransformGwasObject <- function(gwas.obj,marker.selection=1:gwas.obj$N,modify.external=T) { # ,ind.selection=1:gwas.obj$n
#browser()
# OBJECTIVE :
# *
# *
# INPUT :
# * gwas.obj, GWAS.object with 'markers' field being a p x n matrix or data.frame (n=number of genotypes)
# *  modify.external : if TRUE, re-create csv- and bin files of the genotypic data
# *
# OUTPUT :
# *
# *
# *
  if (nrow(gwas.obj$markers)>0) {  # arabidopsis data
    new.gwas.obj <- list()
    new.gwas.obj$description     <- gwas.obj$description
    new.gwas.obj$markers         <- gwas.obj$markers[marker.selection,] # ind.selection
    new.gwas.obj$pheno           <- gwas.obj$pheno
    new.gwas.obj$map             <- gwas.obj$map[marker.selection,]
    #new.gwas.obj$kinship         <- gwas.obj$kinship
    new.gwas.obj$kinship         <- IBS(gwas.obj$markers)
    new.gwas.obj$kinship.asreml  <- gwas.obj$kinship.asreml
    new.gwas.obj$genes           <- gwas.obj$genes
    new.gwas.obj$plant.names     <- gwas.obj$plant.names
    new.gwas.obj$N               <- length(marker.selection)
    new.gwas.obj$n               <- gwas.obj$n # length(ind.selection)
    new.gwas.obj$nchr            <- length(unique(new.gwas.obj$map$chromosome))
    new.gwas.obj$chromosomes     <- sort(unique(new.gwas.obj$map$chromosome))

    if (!is.null(gwas.obj$pc1)) {new.gwas.obj$pc1 <- gwas.obj$pc1}

    #
    new.chr.lengths     <- as.numeric(table(new.gwas.obj$map$chromosome))
    if (new.gwas.obj$nchr > 1) {
      chr.pos         <- c(0,cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)])
      chr.lengths.bp  <- c(0,new.gwas.obj$map$position[cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)]])
      #cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
    } else {
      #cumpositions    <- map$position
      #chr.pos         <- 0
      chr.lengths.bp  <- 0#nrow(marker.object)
    }
    new.gwas.obj$chr.lengths.bp     <- chr.lengths.bp
    #
    if (modify.external) {
      new.gwas.obj$external        <- gwas.obj$external
      new.gwas.obj$external$csv.name <- paste("modified_",gwas.obj$external$csv.name,sep="")
      new.gwas.obj$external$bin.name <- paste("modified_",gwas.obj$external$bin.name,sep="")
      #
      MakeCsv(markers=new.gwas.obj$markers,plant.names=new.gwas.obj$plant.names,file.name=new.gwas.obj$external$csv.name)
      MakeBin(kinship=new.gwas.obj$kinship,plant.names=gwas.obj$plant.names,pheno=gwas.obj$pheno,
              csv.file.name=new.gwas.obj$external$csv.name,bin.file.name=new.gwas.obj$external$bin.name)    # old argument : markers=new.gwas.obj$markers,
      #MakeBin(markers=new.gwas.obj$markers,kinship.file=new.gwas.obj$external$kinship.name,csv.file.name=new.gwas.obj$external$csv.name,bin.file.name=new.gwas.obj$external$bin.name)
      #
    }
    new.gwas.obj$pca             <- gwas.obj$pca

    new.gwas.obj$markerInfo      <- gwas.obj$markerInfo
    new.gwas.obj$real.effects    <- gwas.obj$real.effects
    #new.gwas.obj$     <- gwas.obj$

  } else { # WTCCC-data

    new.gwas.obj <- list()
    #new.gwas.obj$description     <- gwas.obj$description
    new.gwas.obj$pheno           <- gwas.obj$pheno
    new.gwas.obj$map             <- gwas.obj$map[marker.selection,]
    new.gwas.obj$plant.names     <- gwas.obj$plant.names
    new.gwas.obj$N               <- length(marker.selection)
    new.gwas.obj$n               <- gwas.obj$n # length(ind.selection)
    new.gwas.obj$nchr            <- length(unique(new.gwas.obj$map$chromosome))
    new.gwas.obj$chromosomes     <- sort(unique(new.gwas.obj$map$chromosome))
    #
    new.chr.lengths     <- as.numeric(table(new.gwas.obj$map$chromosome))
    if (new.gwas.obj$nchr > 1) {
      chr.pos         <- c(0,cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)])
      chr.lengths.bp  <- c(0,new.gwas.obj$map$position[cumsum(new.chr.lengths)[1:(new.gwas.obj$nchr-1)]])
      #cumpositions    <- map$position + rep(cumsum(chr.lengths.bp),times=chr.lengths)
    } else {
      #cumpositions    <- map$position
      #chr.pos         <- 0
      chr.lengths.bp  <- 0#nrow(marker.object)
    }
    new.gwas.obj$chr.lengths.bp     <- chr.lengths.bp
    #
    new.gwas.obj$external        <- gwas.obj$external
    new.plink.name                  <- paste("modified_",gwas.obj$external$plink.name,sep="")
    new.gwas.obj$external$plink.name <- new.plink.name
    #
    #if (!all(file.exists(paste(new.plink.name,".bim",sep="")),file.exists(paste(new.plink.name,".bed",sep="")),file.exists(paste(new.plink.name,".fam",sep="")))) {
      temp.snp.file <- "tempsnpfile.txt"
      write.table(gwas.obj$map$snp.name[marker.selection],file=temp.snp.file,quote=F,col.names=F,row.names=F)
      system(paste("plink --bfile",gwas.obj$external$plink.name,"--make-bed --extract",temp.snp.file,"--out",new.plink.name))
    #}
    system(paste("plink --bfile",new.plink.name,"--recodeA --out",new.plink.name))
    new.gwas.obj$markers <- t(read.table(file=paste(new.plink.name,".raw",sep=""),header=T)[,-(1:6)])#t(read.table(file="modified_CD_and_wtc1.raw",header=T)[,-(1:6)])
    new.gwas.obj$pca             <- gwas.obj$pca
  }
  if (length(gwas.obj$groups)>0) {
    group.selection <- apply(gwas.obj$groups$position.matrix,1,FUN=function(x){sum(marker.selection %in% (x[1]):(x[2]))})
    G.temp <- gwas.obj$groups$position.matrix[group.selection>0,]
    G.temp[1,1] <- marker.selection[1]
    G.temp[nrow(G.temp),2] <- marker.selection[length(marker.selection)]
    G.temp <- G.temp - G.temp[1,1] + 1
    new.gwas.obj <- AddGroupsToGwasObj(G.temp,new.gwas.obj)
  }
return(new.gwas.obj)
}

