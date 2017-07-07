######################################################################################
########    DRAWING THE FIGURE OF QTLS AS IN Millet et al. 2016        ##############
######################################################################################
##### E. Millet, 27/01/2017

rm(list=ls())
# set working directory
setwd("/home/millet/Documents/PhD-1A/GWAS/Z_FigureGWAS_multiTrait")
# loading package
library(ggplot2)
library(ReporteRs) # R package to format the outputs using Word or Powerpoint
#
######                      ###### 
###### Dataframe formatting ###### 
######                      ###### 

### Dataframe format should be output of the "listQTL" function
# or should at least contain SNP name, snp position, allelic effect, chromosome, trait and/or environment
# MANDATORY : 
### column containing chromosome should be named "Chromosome" 
### column containing trait should be named "Trait" (either phenotypic trait or environment)
### column containing SNP effect should be named "effect" 
### column containing SNP position should be named "Marker_Position"
# Example 1 : different phenotypic trait
listeSNP <- read.table("ListeQTL_LOD45_DROPS_Response.csv",sep=",",he=TRUE)
# Example 2: different environment
listeSNP <- read.table("Fixed_modelMEML_33+15QTL_GY_perenvt_effB73.csv",sep=",",he=TRUE)
# center and reduce the allelic effect (because of the different units)
listeSNP$eff <- sapply(1:nrow(listeSNP), function(x) 
  (listeSNP$effect[x]-mean(listeSNP$effect[listeSNP$Trait==listeSNP$Trait[x]],na.rm=T))/sd(listeSNP$effect[listeSNP$Trait==listeSNP$Trait[x]],na.rm=T) )
listeSNP$eff <- listeSNP$effect

# OPTION : adding a sorting column, to order the trait or environment on the graph
# Example 1 : no sorting
listeSNP$tri <- 1
# Example 2 : sorting
tri <- read.table("TriClusterMaxTemp.csv",sep=",",he=T)
listeSNP$tri<-tri$triscenariowatertemp[match(listeSNP$Trait,tri$Envtreduce)]


### Adding lines to the dataframe that contain the chromosome limits
# otherwise the chromosome length on the figure will depend on the QTLs only
# (i.e. if a chromosome display no QTL it would not be drawn)
limits <- read.table("/home/millet/Documents/PhD-1A/GWAS/A_50K+600K_GyGnbGsAnth/ChrLimits_RefGenV2.csv",sep=",",he=T)
#
df           <- data.frame(matrix(ncol=ncol(listeSNP),nrow=20 ))           # empty datafram with 20 lines (2 per chromosome) and as many columns as the QTL dataframe
names(df)    <- names(listeSNP)                                            # columns names = names of the QTL dataframe
df$Trait     <- levels(factor(listeSNP$Trait))[1]                          # give 1 trait value (any one)
df$tri       <- listeSNP$tri[listeSNP$Trait==unique(df$Trait)][1]          # (option) the order number of the trait or situation
df[,c("Chromosome","Marker_Position")] <- limits[,c(1,2)]                  # add the physical limits of the chromosomes
#
dat          <- rbind(listeSNP,df)                                         # concatenate this small dataframe with the QTL one

### Adding a column with the allelic effect direction (for points color)
dat$coul <- ifelse(dat$eff>0,"pos","neg")
dat <- droplevels(dat)


######        ###### 
######  Plot  ###### 
######        ###### 

### Loading the bin position to add vertical lines to the figure
### WARNINGS : particular case of the maize genome ! 
### other species do not have bins
dat.vline<-read.table("/home/millet/Documents/PhD-1A/GWAS/A_50K+600K_GyGnbGsAnth/PosBinRefGenV2.csv",sep=",",he=T)
names(dat.vline)[1] <- "Chromosome"

### Loading the ggplot theme ready for the figure
# this one is handmade and fit nicely the figure
# for more information check http://docs.ggplot2.org/current/theme.html
source("/home/millet/Documents/PhD-1A/ScriptR/themeggplot_GWAS.R")

### Draw the plot

# labels of the yaxis
y.labels <- "Traits"       # example 1
y.labels <- "Environments" # example 2

qtlplot <-                                                   # R object containing the plot
  ggplot(dat,                                                # data set
         aes(x      = Marker_Position,                       # x data
             y      = reorder(Trait,-tri),                   # y data sorted as 'tri' ('reorder' function)
             size   = factor(abs(eff)),                      # points size proportionnal to allelic effect (absolute value)
             colour = factor(coul) ) ) +                     # points color depends on the effect direction
         my_theme_empty_GWAS +                               # use the 'handmade' ggplot theme 
         ylab(y.labels) + xlab("Chromosomes")  +             # labelling the plot axis
         geom_vline(aes(xintercept = position),              # add vertical lines at the bin position
                   data = dat.vline,                         # the bin positions are contained in the 'dat.vline' dataframe
                   linetype=1,                               # type of line (1 = full)
                   color = "white")   +                      # color of line
        geom_point(alpha = I(0.7))  +                        # add the points with a slight transparency in case of overlap
        facet_grid(". ~ Chromosome",                         # split of the plot according to the chromosomes on the x axis
                    scales = "free",                         # do not resize the x axise (otherwise every chromosome has the same size)
                    space = "free",                          # do not add extra space between facet
                   switch="both") +                          # place the chromosome labels at the bottom
        scale_colour_manual("coul",                          # manually ascribe a color to the allelic direction (column 'coul')...  
                            labels = c("neg", "pos"),        # ... which could be either "neg" or "pos" ...
                            values = c("darkblue","green4")) # ... and that has a corresponding color 
  
  
######                                      ###### 
######  Saving the figure using in .pptx    ###### 
######                                      ###### 

### Loading an empty powerpoint template 
# this one is a .pptx file, containing one empty slide
# for more information check http://davidgohel.github.io/ReporteRs/
mydoc = pptx(template = 'templatevide_ReporteRs.pptx')  
mydoc = addSlide( mydoc, "Titre et contenu" )     # adding a new slide (necessary despite the first empty slide)
mydoc = addPlot( doc = mydoc,                     # adding a plot the the document
                 fun = print,                     # "print" the graph on the slide 
                 x = qtlplot,                     # x = R object containing the plot
                 offx = 1, offy = 1,              # position of the plot on the slide
                 width = 10, height = 8 )         # plot size
mydoc = addDate(mydoc)                            # adding date on the slide
writeDoc( doc=mydoc, file="GraphQTL.pptx")        # writing the .pptx file
