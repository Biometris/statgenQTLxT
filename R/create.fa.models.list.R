create.fa.models.list <- function(m.G,m.E,het.G,het.E) {
  d <- data.frame(m.G=m.G,m.E=m.E,het.G=het.G,het.E=het.E)
  rownames(d) <- paste0('model',1:nrow(d))
return(d)
}
