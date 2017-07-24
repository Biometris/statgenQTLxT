######### to be tested !!

#### clarify the role of B, checks ....

BLUP_function <- function(Y,K,X=matrix(rep(1,nrow(K))),B=NULL,Cm,Dm) {

# X default used to be : X=data.frame()

# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed

  stopifnot(nrow(Y)==nrow(K))
  stopifnot(ncol(K)==nrow(K))

  nc    <- ncol(X)

  if (nc > 0) {
    stopifnot(nrow(X)==nrow(K))
    stopifnot(!is.null(B))
    stopifnot(class(B) %in% c('matrix','data.frame'))
    stopifnot(nrow(B)==nc)
  }

  n <- ncol(K)

  p <- ncol(Y)

  R <- solve(K)

  w <- eigen(R)

  Lambda.R <- diag(w$values)

  U <- w$vectors

  ###

  Dm.sqrt.inv <- MatrixRoot(ginv(Dm))

  w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

  Q.1 <- w$vectors

  Lambda.1 <- w$values

  if (nc > 0) {

    S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% (Y - X %*% B) %*% MatrixRoot(Dm) %*% Q.1)

  } else {

    S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% Y %*% MatrixRoot(Dm) %*% Q.1)

  }

  Dm.sqrt.inv <- MatrixRoot(solve(Dm))

  w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

  Q.1 <- w$vectors

  mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)

return(mu)
}
