cross.validate.FA.models <- function(Y,K,X,
                                     fa.models,fold.obj,
                                     max.diag=1000,
                                     max.iter.em=1000,tol.em=0.0001,
                                     max.iter.em.cv=1000,tol.em.cv=0.0001) {

#max.diag=1000;max.iter.em=1000;tol.em=0.0001;max.iter.em.cv=1000;tol.em.cv=0.0001;m=1;j=1

  n.fold   <- length(fold.obj)
  n.models <- nrow(fa.models)

  cv.results.main           <- matrix(NA,n.models,ncol(Y))
  rownames(cv.results.main) <- rownames(fa.models)
  colnames(cv.results.main) <- colnames(Y)

  LL <- rep(NA,n.models)

  for (m in 1:n.models) {

    full.set <- EM_function_FA(Y=Y,K=K,X=X,
                      max.iter.em=max.iter.em,
                      tol.em=tol.em,
                      m.G=fa.models$m.G[m],
                      m.E=fa.models$m.E[m],
                      max.diag=max.diag,
                      Cm.het=fa.models$het.G[m],
                      Dm.het=fa.models$het.E[m])

    LL.m <- full.set$log.lik
    LL[m] <- LL.m

    cv.results <- matrix(NA,n.fold,ncol(Y))
    rownames(cv.results) <- paste0('set',1:n.fold)
    colnames(cv.results) <- colnames(Y)

    for (j in 1:n.fold) {
      cv.j <- EM_function_FA(Y=Y[fold.obj[[j]]$train,],
                             K=K[fold.obj[[j]]$train,fold.obj[[j]]$train],
                             X=as.matrix(X[fold.obj[[j]]$train,]),
                             max.iter.em=max.iter.em.cv,
                             tol.em=tol.em.cv,
                             m.G=fa.models$m.G[m],
                             m.E=fa.models$m.E[m],
                             max.diag=max.diag,
                             Cm.het=fa.models$het.G[m],
                             Dm.het=fa.models$het.E[m],
                             Cm.start=full.set$Cm,
                             Dm.start=full.set$Dm)

      # G-BLUP for the test set
      Gpred <- t(matrix(cv.j$pred$predicted,ncol=cv.j$n) %*% solve(K[fold.obj[[j]]$train,fold.obj[[j]]$train]) %*% K[fold.obj[[j]]$train,fold.obj[[j]]$test])
      colnames(Gpred) <- colnames(Y)
      rownames(Gpred) <- rownames(Y)[fold.obj[[j]]$test]

      # differences
      #Y[rownames(Gpred),colnames(Gpred)] - Gpred

      # correlations, per trait
      cv.j.cor <- mapply(cor,as.data.frame(Y[rownames(Gpred),colnames(Gpred)]),as.data.frame(Gpred))
      cv.j.cor[cv.j.cor <0 ] <- 0
      cv.results[j,]         <- cv.j.cor
    }
    cv.results.mean     <- apply(cv.results,2,mean)
    cv.results.main[m,] <- cv.results.mean
  }

return(data.frame(cv.results.main,LL=LL))
}
