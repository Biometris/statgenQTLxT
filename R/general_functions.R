







##########










# matrix.to.par(matrix(1:9,ncol=3))











#x <- as.numeric(t(matrix(as.numeric(GWAS.obj$markers[mrk,rownames(Y)]))))
#newX <- rbind(X,t(matrix(x)) %*% Uk) #t(as.matrix(data.frame(intercept=rep(1,nrow(Y)),marker=x)))
#V.inv.array <- make.V.inv.array(Vg=Vg,Ve=Ve,Dk=Dk)


#####################################################
















##############################################################################################################

    #LL0 <- LL.diag(Y=Y-fitted.mean0,V.array=V.array,V.inv.array=V.inv.array)

  ##################

  #LL1 <- LL.diag(Y=Y-fitted.mean1,V.array=V.array,V.inv.array=V.inv.array)
  #LRT <- 2 * (LL1 - LL0)

  ##############

  #LRT.pvalue <- pchisq(LRT,df=df.red-df.full,lower.tail=F)

  #result <- c(LRT.pvalue,LRT,LL1,LL0)
  #names(result) <- c('pvalue','LRT','LL1','LL0')
  #return(result)

  # est1$effects.estimates[-(1:p)]
  # est1$effects.sd[-(1:p)]

  #  return(list(pvalue=F.pvalue,LRT=LRT,LL1=LL1,LL0=LL0,df=df.red-df.full,
