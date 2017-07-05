bic.FA.models <- function(Y,K,X,
                          fa.models,
                          max.diag=100,
                          max.iter.em=1000,
                          tol.em=0.0001,
                          stopifdecreasing=FALSE) {

#max.diag=1000;max.iter.em=1000;tol.em=0.0001;max.iter.em.cv=1000;tol.em.cv=0.0001;m=1;j=1

  n.models <- nrow(fa.models)

  LL.E <- LL <- n.iter <- conv <- decreased <- rep(NA,n.models)

  for (m in 1:n.models) {

    cat('Model', m,'\n')

    full.set <- EM_function_FA(Y=Y,K=K,X=X,
                      max.iter.em=max.iter.em,
                      tol.em=tol.em,
                      m.G=fa.models$m.G[m],
                      m.E=fa.models$m.E[m],
                      max.diag=max.diag,
                      Cm.het=fa.models$het.G[m],
                      Dm.het=fa.models$het.E[m],
                      compute.log.lik=TRUE,
                      stopifdecreasing=stopifdecreasing)

    LL[m]        <- full.set$log.lik2
    LL.E[m]      <- full.set$log.lik
    n.iter[m]    <- full.set$n.iter
    conv[m]      <- full.set$converged
    decreased[m] <- full.set$decreased
  }

  fa.models$het.G <- as.numeric(fa.models$het.G)
  fa.models$het.E <- as.numeric(fa.models$het.E)

  p <- ncol(Y)
  n <- nrow(Y)
  nc<- ncol(X)

  npar <- (as.numeric(apply(fa.models,1,sum)) + nc) * p

  bic <- -2 * LL + log(n * p) * npar

  aic <- -2 * LL + 2 * (n*p / (n*p - npar - 1)) * npar

return(data.frame(LL=LL,LL.E=LL.E,bic=bic,aic=aic,n.iter=n.iter,converged=conv,decreased=decreased))
}
