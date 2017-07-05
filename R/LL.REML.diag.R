LL.REML.diag <- function(Y,X,V.array,V.inv.array) {

# X is the c x n covariate matrix, c being the number of covariates and
#      n being the number of genotypes
#      c has to be at least one (typically an intercept)
# Y is the p x n matrix of observed phenotypes, on p traits or environments

#
# No missing values are allowed in any of X,Y and x
#
# It is assumed that X,Y and x have already been rotated by Uk, where Uk is such that
# the kinship matrix equals K = Uk %*% Dk %*% t(Uk)
# (R: w <- eigen(K); Dk <- diag(w$values); Uk <- w$vectors)
# Next, the original X, Y and x are post-multiplied by Uk, e.g. Y <- Y %*% Uk;
# see Zhou and Stephens 2014, supplement.
# It is the rotated versions that are the input of this function.

  stopifnot(ncol(Y)==ncol(X))

  nc      <- nrow(X)
  y       <- matrix(Y)
  p       <- nrow(Y)
  n       <- ncol(Y)

  q.scal <- sum(sapply(1:n,function(i){t(matrix(Y[,i])) %*% V.inv.array[i,,] %*% matrix(Y[,i])}))

  if (nc > 0) {
    q.vec  <- matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]),V.inv.array[i,,] %*% matrix(Y[,i]))})),1,sum))
    Q.matrix <-  matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]) %*% t(matrix(X[,i])),V.inv.array[i,,])})),1,sum),ncol=p)
  } else {
    q.vec  <- matrix(rep(0,p))
    Q.matrix <- matrix(rep(0,p*p),ncol=p)
  }

  det.part <-  -0.5 * (sum(sapply(1:n,function(i){log(det(V.array[i,,]))})))
  det.part2 <- 0.5 * p * log(det(X %*% t(X)))
  det.part3 <- -0.5 * log(det(Q.matrix))

  quad.form.part <-  -0.5 * (q.scal - as.numeric(t(q.vec) %*% solve(Q.matrix) %*% q.vec))

  REML.log.lik <- quad.form.part + det.part + det.part2 + det.part3

return(REML.log.lik)
}
# LL.REML.diag(Y=Yt,X=X,V.array=V.array,V.inv.array=V.inv.array)

