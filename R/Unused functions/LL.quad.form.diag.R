#######################################################################################################
# the same as LL.diag, but without the det-part

LL.quad.form.diag <- function(Y,X=data.frame(),V.inv.array) {

# X is the c x n covariate matrix, c being the number of covariates and
#      n being the number of genotypes
#      c has to be at least one (typically an intercept)
# Y is the p x n matrix of observed phenotypes, on p traits or environments
#
# No missing values are allowed in any of X,Y and x
#
# It is assumed that X,Y have already been rotated by Uk, where Uk is such that
# the kinship matrix equals K = Uk %*% Dk %*% t(Uk)
# (R: w <- eigen(K); Dk <- diag(w$values); Uk <- w$vectors)
# Next, the original X, Y are post-multiplied by Uk, e.g. Y <- Y %*% Uk;
# see Zhou and Stephens 2014, supplement.
# It is the rotated versions that are the input of this function.

  #stopifnot(ncol(Y)==ncol(X))

  # Y=Yt;X=data.frame()
  # Y=Yt;X=rbind(X,x)
  # Y=Yt-fitted.mean # X=data.frame()

  # Y=Yt;X=Xt;V.inv.array=V.inv.array

  nc      <- nrow(X)
  y       <- matrix(Y)
  p       <- nrow(Y)
  n       <- ncol(Y)

  q.scal <- sum(sapply(1:n,function(i){t(matrix(Y[,i])) %*% V.inv.array[i,,] %*% matrix(Y[,i])}))

  quad.form.part <-  q.scal

  if (nc > 0) {

    if (p==1 & nc==1) {

      q.vec  <- matrix(sum(sapply(1:n,function(i){kronecker(matrix(X[,i]),V.inv.array[i,,] %*% matrix(Y[,i]))})))
      Q.matrix <-  matrix(sum(sapply(1:n,function(i){kronecker(matrix(X[,i]) %*% t(matrix(X[,i])),V.inv.array[i,,])})))
      quad.form.part <-  quad.form.part -1 * (as.numeric(t(q.vec) %*% solve(Q.matrix) %*% q.vec))

    } else {

      q.vec  <- matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]),V.inv.array[i,,] %*% matrix(Y[,i]))})),1,sum))
      Q.matrix <-  matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]) %*% t(matrix(X[,i])),V.inv.array[i,,])})),1,sum),ncol=p*nc)
      quad.form.part <-  quad.form.part -1 * (as.numeric(t(q.vec) %*% solve(Q.matrix) %*% q.vec))
    }

  }


return(quad.form.part)
}
# LL.diag(Y=Yt,X=X,V.array=V.array,V.inv.array=V.inv.array)
