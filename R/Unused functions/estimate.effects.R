estimate.effects <- function(X,x=numeric(),Y,Dk,V.inv.array,return.all.effects=T) {
# return.all.effects : if TRUE, also return the effect estimates and standard errors
# for the other covariates (i.e. those in X)

# X is the c x n covariate matrix, c being the number of covariates and
#      n being the number of genotypes
#      c has to be at least one (typically an intercept)
# Y is the p x n matrix of observed phenotypes, on p traits or environments
# x (optional) is an additional covariate.
#
# No missing values are allowed in any of X,Y and x
#
# It is assumed that X,Y and x have already been rotated by Uk, where Uk is such that
# the kinship matrix equals K = Uk %*% Dk %*% t(Uk)
# (R: w <- eigen(K); Dk <- diag(w$values); Uk <- w$vectors)
# Next, the original X, Y and x are post-multiplied by Uk, e.g. Y <- Y %*% Uk;
# see Zhou and Stephens 2014, supplement.
# It is the rotated versions that are the input of this function.
#
# The matrix Uk is not needed as input, but the diagonal matrix Dk is.
# We also need the n x p x p dimensional array V.inv.array, which can be obtained
# from the function make.V.inv.array (which in turn requires Dk and the current values of Vg and Ve).
# V.inv.array contains for each genotype i=1,...,n the p x p matrix V_{l}^{-1}
# (in the notation of Zhou and Stephens)
#
# If return.all.effects==T, only the p effect estimates for the additional covariate are returned,
# (the first p * c are actually computed, but not returned)
# If no x is provided, return.all.effects is automatically put to TRUE

  # stopifnot(...)
  # check missing values
  # comment on missing values
  # to do : partial update...

# X=Xt;Y=Yt;Dk=Dk;V.inv.array=V.inv.array;return.all.effects=T ; x=numeric()
# X=Xt;Y=Yt;Dk=Dk;V.inv.array=V.inv.array;return.all.effects=T ; x=xt
# X=X;Y=Yt;Dk=Dk;V.inv.array=V.inv.array;return.all.effects=T ; x=x
# return.all.effects=T

  # the existing ... is already assumed to be post-multiplied by Uk
  nc   <- nrow(X)
  p    <- nrow(Y)

  if (length(x)==0) {
    return.all.effects <- T
    nc.total <- nc
  } else {
    nc.total <- nc + 1
  }

  if (length(x) > 0) {X    <- rbind(X,t(matrix(as.numeric(x))))}
  # the last p coefficients should be the marker-effects,
  #  the first p*nc should correspond to the other coefficients

  if (p==1 & nc.total==1) {
    Vbeta <- matrix(sum(sapply(1:n,function(m){kronecker((matrix(X[,m])) %*% t(matrix(X[,m])),V.inv.array[m,,])})))
  } else {
    Vbeta <- matrix(apply(sapply(1:n,function(m){kronecker((matrix(X[,m])) %*% t(matrix(X[,m])),V.inv.array[m,,])}),1,sum),ncol=p*(nc.total))
  }

  M     <- solve(Vbeta)

  M.sub <- NULL
  wald  <- NA

  if (length(x) > 0) {
    M.sub <- solve(as.matrix(Vbeta[-(1:(p*nc)),-(1:(p*nc))]))
  }

  if (p==1 & nc.total==1) {
    v     <- sum(sapply(1:n,function(m){kronecker((matrix(X[,m])),V.inv.array[m,,] %*% matrix(Yt[,m]))}))
  } else {
    v     <- apply(sapply(1:n,function(m){kronecker((matrix(X[,m])),V.inv.array[m,,] %*% matrix(Yt[,m]))}),1,sum)
  }

  v <- matrix(v)

  if (return.all.effects) {
    effects.estimates <- as.numeric(M %*% v)
    effects.sd        <- sqrt(diag(M))
    #wald              <- as.numeric(t(matrix(effects.estimates[-(1:(p*nc))])) %*% M[-(1:(p*nc)),-(1:(p*nc))] %*% matrix(effects.estimates[-(1:(p*nc))]))
    if (length(x) > 0) {
      wald              <- as.numeric(t(matrix(effects.estimates[-(1:(p*nc))])) %*% M.sub %*% matrix(effects.estimates[-(1:(p*nc))]))
    }
  } else {
    effects.estimates <- as.numeric(M %*% v)[-(1:(p*nc))]
    effects.sd        <- sqrt(diag(M))[-(1:(p*nc))]
    #wald              <- as.numeric(t(matrix(effects.estimates)) %*% M[-(1:(p*nc)),-(1:(p*nc))] %*% matrix(effects.estimates))
    if (length(x) > 0) {
      wald              <- as.numeric(t(matrix(effects.estimates)) %*% M.sub %*% matrix(effects.estimates))
    }
  }

return(list(effects.estimates=effects.estimates,effects.sd=effects.sd,wald=wald))
}
#
