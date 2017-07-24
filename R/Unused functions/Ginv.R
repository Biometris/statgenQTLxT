Ginv <- function(M) {

  Minv <- try(MASS::ginv(M), silent=TRUE)

  if (class(Minv)=="try-error") {
    Minv <- solve(M)
  }

return(Minv)
}
