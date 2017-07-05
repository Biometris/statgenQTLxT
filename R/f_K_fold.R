# taken from :
# http://stackoverflow.com/questions/7402313/generate-sets-for-cross-validation-in-r
# Code by wojciech sobala
f_K_fold <- function(Nobs, K=5){
    id <- sample(rep(seq.int(K), length.out=Nobs))
    l <- lapply(seq.int(K), function(x) list(
         train = which(x!=id),
         test  = which(x==id)
    ))
    return(l)
}