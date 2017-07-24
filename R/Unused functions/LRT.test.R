LRT.test <- function(Y,X,x,Dk,V.inv.array,SS0=NULL) { # V.array=V.array
  # Y=Yt;X=Xt;x=xt
  # Y=Yt;X=Xt;x=xt;Dk=Dk;V.inv.array=V.inv.array;SS0=SS0
  # Y=Yt;X=Xt.red;x=xt;Dk=Dk;V.inv.array=V.inv.array.red;SS0=SS0.red
  nc <- nrow(X)

  n <- ncol(X)

  p <- nrow(Y)

  df.full <- (n - (nc+1))*p

  df.red  <- (n - nc)*p

  ### Null model with the trait specific means only
  #    (which should be the model for which the Vg and Ve estimates were obtained)

  if (is.null(SS0)) {
    est0 <- estimate.effects(X=X,Y=Y,Dk=Dk,V.inv.array=V.inv.array,return.all.effects=T)

    fitted.mean0 <- matrix(est0$effects.estimates,ncol=length(est0$effects.estimates)/p) %*% X

    SS0 <- LL.quad.form.diag(Y=Y-fitted.mean0,V.inv.array=V.inv.array)
  }
  ### alternative model, with additional trait specific marker effect; calculate for each marker

  est1 <- estimate.effects(x=x,X=X,Y=Y,Dk=Dk,V.inv.array=V.inv.array,return.all.effects=T)

  fitted.mean1 <- matrix(est1$effects.estimates,ncol=length(est1$effects.estimates)/p) %*% rbind(X,x)

  SS1 <- LL.quad.form.diag(Y=Y-fitted.mean1,V.inv.array=V.inv.array)

  Fstat <- ((SS0-SS1) / (SS1)) * (df.full) / (df.red-df.full)

  F.pvalue <- pf(q=Fstat, df1=df.red-df.full, df2=df.full, lower.tail = F)

  return(list(pvalue=F.pvalue,Fstat=Fstat,
              SS1=SS1,SS0=SS0,
              df=df.red-df.full,
              effects=est1$effects.estimates[-(1:(nc*p))],
              effects.se=est1$effects.sd[-(1:(nc*p))],
              wald=est1$wald
              ))
}
