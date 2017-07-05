Ginv <- function(M) {

  Minv <- try(ginv(M), silent=TRUE)

  if (class(Minv)=="try-error") {
    Minv <- solve(M)
  }

return(Minv)
}