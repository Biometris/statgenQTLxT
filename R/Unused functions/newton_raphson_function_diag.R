newton_raphson_function_diag <- function(Y,K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,
                        tol.nr=0.0001,max.iter.nr=300,REML=F,
                        Cm.start=NULL,Dm.start=NULL,nr.scaling=1) {

# Cm.start=test.em$Cm;Dm.start=test.em$Dm;X=matrix(rep(1,nrow(K)));Dm.is.diagonal=F
# X=matrix(rep(1,nrow(K))); REML=TRUE;Cm.start=test.em$Cm;Dm.start=test.em$Dm; nr.scaling=1; Dm.is.diagonal=F


# to do : allow for Dm.is.diagonal=TRUE

# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                with row and column names,
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed

# Dm.is.diagonal:

# tol.nr       : tolerance in the newton-raphson-algorithm : stop when the increase in log-lik. is less than tol
# max.iter.enr  : maximum number of iterations in the newton-raphson algorithm
# Cm.start     : starting values for Cm=solve(Vg), as p x p matrix
# Dm.start     : starting values for Dm=solve(Ve), as p x p matrix

  # also check : missing data

  # X=matrix(rep(1,nrow(K)))
  stopifnot(nrow(Y)==nrow(K))
  stopifnot(ncol(K)==nrow(K))

  nc    <- ncol(X)
  if (nc > 0) {stopifnot(nrow(X)==nrow(K))}

  n <- ncol(K)

  p <- ncol(Y)

  ###############################################################

  if (is.null(Cm.start)) {
    Cm <- ginv((cor(Y) + diag(p)) / 4) #0.5 * cor(Y)
  } else {
    Cm <- Cm.start
  }

  if (is.null(Dm.start)) {
    Dm <- diag(p) / 2
  } else {
    Dm <- Dm.start
  }

  par.vector <- rbind(matrix.to.par(solve(Cm)),matrix.to.par(solve(Dm)))

  ###############################################################

  w <- eigen(K)

  Dk <- diag(w$values)

  Uk <- w$vectors  #sum(abs(K - Uk %*% Dk %*% t(Uk)))

  rownames.Y <- rownames(Y)

  Yt <- t(Y) %*% Uk

  colnames(Yt) <- rownames.Y

  if (nc > 0) {
    X <- t(X) %*% Uk
  }

  continue <- T

  iter <- 1

  # Q is the number of parameters in each of Cm and Dm (and their inverses).
  # The total number of parameters (apart from the fixed effects) is 2*Q
  Q <- p * (p + 1) / 2

  # starting value for the log-lik of the previous iteration
  NR.log.lik.old <- 10000000000

  while (continue & iter < max.iter.nr) {

    # check if the current parameter values imply positive definite Vg, Ve,
    # and if not make them pos. definite
    if (!is.positive.definite(par.to.matrix(par.vector[1:Q,]))) {
      par.vector[1:Q,] <- matrix.to.par(nearPD(par.to.matrix(par.vector[1:Q,]),keepDiag=F)$mat)
    }
    if (!is.positive.definite(par.to.matrix(par.vector[Q + 1:Q,]))) {
      par.vector[Q + 1:Q,] <- matrix.to.par(nearPD(par.to.matrix(par.vector[Q +  1:Q,]),keepDiag=F)$mat)
    }

    V.inv.array  <- make.V.inv.array(Vg=par.to.matrix(par.vector[1:Q,]),
                                     Ve=par.to.matrix(par.vector[Q + 1:Q,]),
                                     Dk=Dk)

    V.array      <- make.V.array(Vg=par.to.matrix(par.vector[1:Q,]),
                                 Ve=par.to.matrix(par.vector[Q + 1:Q,]),
                                 Dk=Dk)

    if (REML) {
      NR.log.lik <- LL.REML.diag(Y=Yt,X=X,V.array=V.array,V.inv.array=V.inv.array)
    } else {
      NR.log.lik <- LL.diag(Y=Yt,X=X,V.array=V.array,V.inv.array=V.inv.array)
    }

    if (abs(NR.log.lik - NR.log.lik.old) < tol.nr) {continue <- F}

    if (REML) {

      GR <- LL.REML.grad.diag(Y=Yt,X=X,Dk=Dk,V.inv.array=V.inv.array)

      grad       <- GR$grad

      Q.matrix   <- GR$Q.matrix

      q.vec      <- GR$q.vec

      Q.G.output <- GR$Q.G.output

      Q.E.output <- GR$Q.E.output

      Q.g.output <- GR$Q.g.output

      Q.e.output <- GR$Q.e.output

    } else {

      GR <- LL.grad.diag(Y=Yt,X=X,Dk=Dk,V.inv.array=V.inv.array)

      grad       <- GR$grad

      Q.matrix   <- GR$Q.matrix

      q.vec      <- GR$q.vec

      Q.G.output <- GR$Q.G.output

      Q.E.output <- GR$Q.E.output

      Q.g.output <- GR$Q.g.output

      Q.e.output <- GR$Q.e.output

    }

    if (REML) {
      hess <- LL.REML.hess.diag(Y=Yt,X=X,Dk=Dk,V.inv.array=V.inv.array,
                           Q.G.output=Q.G.output,Q.E.output=Q.E.output,
                           Q.g.output=Q.g.output,Q.e.output=Q.e.output,
                           q.vec=q.vec,Q.matrix=Q.matrix
                           )
    } else {
      hess <- LL.hess.diag(Y=Yt,X=X,Dk=Dk,V.inv.array=V.inv.array,
                           Q.G.output=Q.G.output,Q.E.output=Q.E.output,
                           Q.g.output=Q.g.output,Q.e.output=Q.e.output,
                           q.vec=q.vec,Q.matrix=Q.matrix
                           )

    }

    par.vector <- par.vector - nr.scaling * ginv(hess) %*% grad

    NR.log.lik.old <- NR.log.lik

    cat('Iteration ',iter,'\n loglik :',NR.log.lik,'\n','parameter-estimates: ',round(as.numeric(par.vector),4),'\n')

    iter <- iter + 1

  }

  Cm <- solve(par.to.matrix(par.vector[1:Q,]))
  Dm <- solve(par.to.matrix(par.vector[Q + 1:Q,]))

return(list(Cm=Cm,Dm=Dm,log.lik=NR.log.lik,tol.nr=tol.nr,converged=(!continue)))

}
