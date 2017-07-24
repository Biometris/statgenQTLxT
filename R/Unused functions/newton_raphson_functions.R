# requires :   source('newton_raphson_functions.R'); general functions

newton_raphson_function <- function(Y,K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,
                        tol.nr=0.0001,max.iter.nr=300,REML=F,
                        Cm.start=NULL,Dm.start=NULL,nr.scaling=1) {

# max.iter.nr <- 5; Cm.start=test.em$Cm;Dm.start=test.em$Dm;X=matrix(rep(1,nrow(K)));Dm.is.diagonal =F; nr.scaling=1


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

  if (nc > 0) {bigX <- kronecker(t(X),diag(p))}

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

    H <- kronecker(Dk,par.to.matrix(par.vector[1:Q,])) + kronecker(diag(n),par.to.matrix(par.vector[Q + 1:Q,]))

    Hinv <- solve(H)

    if (nc > 0) {
      P <- Hinv - Hinv %*% bigX %*% solve(t(bigX) %*% Hinv %*% bigX) %*% t(bigX) %*% Hinv
    } else {
      P <- Hinv
    }

    if (REML) {
      NR.log.lik <- LL.REML(Y=Yt,H=H,P=P,X=X)
    } else {
      NR.log.lik <- LL(Y=Yt,H=H,P=P)
    }

    if (abs(NR.log.lik - NR.log.lik.old) < tol.nr) {continue <- F}

    if (REML) {
      grad <- LL.REML.grad(Y=Yt,P=P,Dk=Dk)
      #system.time(grad <- LL.REML.grad(Y=Yt,P=P,Dk=Dk))  # 2.84    0.07    2.90
      #system.time(grad <- LL.REML.grad.diag(Y=Yt,P=P,Dk=Dk))
    } else {
      grad <- LL.grad(Y=Yt,P=P,Dk=Dk,Hinv=Hinv)
    }

    if (REML) {
      hess <- LL.REML.hess(Y=Yt,P=P,Dk=Dk)
    } else {
      hess <- LL.hess(Y=Yt,P=P,Dk=Dk,Hinv=Hinv)
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




###################################################################################################

# to be completed and corrected
LL.full <- function(Y,X,Cm,Dm,K) {
  n <- ncol(Y)
  p <- nrow(Y)
  stopifnot(p==ncol(Cm) & p==ncol(Dm))
  H <- kronecker(solve(Cm),K) + kronecker(solve(Dm),diag(n))
  Hinv <- solve(H)
  P <- Hinv - Hinv %*% X %*% solve(t(X) %*% Hinv %*% X) %*% t(X) %*% Hinv
}

# The log-likelihood, as function of the trait-matrix Y and H and P.
# Parametrization as in Zhou and Stephens, supplement p.19, equation (26);
# H and P should be computed before as in equations (28-29).
# Given p traits measured on n genotypes, Y needs to be a p x n matrix,
# and H and P need to be np x np
LL <- function(Y,H,P) {
  stopifnot(ncol(Y)==ncol(H)/nrow(Y))
  stopifnot(ncol(Y)==ncol(P)/nrow(Y))
  y <- matrix(Y)
  # optional: include the constant -0.5 * ncol(Y) * nrow(Y) * log(2 * pi)
return(-0.5 * log(det(H)) - 0.5 * t(y) %*% P %*% y)
}

# The REML log-likelihood, as function of the trait-matrix Y and H, P and X.
# Parametrization as in Zhou and Stephens, supplement p.19, equation (27);
# H and P should be computed before as in equations (28-29).
# Given p traits measured on n genotypes, Y needs to be a p x n matrix,
# and H and P need to be np x np
# X needs to be n x c
LL.REML <- function(Y,H,P,X) {
  stopifnot(ncol(Y)==ncol(H)/nrow(Y))
  stopifnot(ncol(Y)==ncol(P)/nrow(Y))
  stopifnot(ncol(Y)==ncol(X))
  nc      <- ncol(X)
  y       <- matrix(Y)
  Hinv    <- solve(H)
  result  <- -0.5 * log(det(H)) - 0.5 * t(y) %*% P %*% y
  logdet1 <- log(det(X %*% t(X)))
  logdet2 <- log(det(kronecker(X,diag(nrow(Y))) %*% Hinv %*% kronecker(t(X),diag(nrow(Y)))))
  result  <- result + 0.5 * logdet1 -0.5 * logdet2
  # optional: include the constant -0.5 * ncol(Y) * nrow(Y) * log(2 * pi)
#return(-0.5 * log(det(H)) - 0.5 * t(y) %*% P %*% y)
return(result)
}
# LL.REML(Y=Yt,H=H,P=P,X=X)


# The gradient of the log-likelihood, as function of the trait-matrix Y,
# the vcov matrix H, projection matrix P, and diagonal matrix Dk, containing the
# eigenvalues of the kinship matrix.
# Parametrization as in Zhou and Stephens, supplement p.19, equation (32-33);
# H and P should be computed before as in equations (28-29).
# Given p traits measured on n genotypes, Y needs to be a p x n matrix,
# and H and P need to be np x np

LL.grad <- function(Y,P,Dk,Hinv) {

  p <- nrow(Y)
  y <- matrix(Y)

  grad.G <- matrix(0,p,p)
  grad.E <- matrix(0,p,p)

  n <- ncol(Dk)

  for (j in 1:p) {
    for (i in 1:j) {
      ij.matrix        <- make.ij.matrix(i=i,j=j,p=p)
      grad.G[i,j]      <-  -0.5 * sum(diag(Hinv %*% kronecker(Dk,ij.matrix)))      + 0.5 * t(y) %*% P %*% kronecker(Dk,     ij.matrix) %*% P %*% y
      grad.E[i,j]      <-  -0.5 * sum(diag(Hinv %*% kronecker(diag(n),ij.matrix))) + 0.5 * t(y) %*% P %*% kronecker(diag(n),ij.matrix) %*% P %*% y
    }
  }

  diag(grad.G) <- diag(grad.G) / 2
  diag(grad.E) <- diag(grad.E) / 2

return(matrix(c(grad.G[upper.tri(grad.G,diag=T)],grad.E[upper.tri(grad.E,diag=T)])))
}

# The gradient of the REML log-likelihood, as function of the trait-matrix Y,
# the vcov matrix H, projection matrix P, and diagonal matrix Dk, containing the
# eigenvalues of the kinship matrix.
# Parametrization as in Zhou and Stephens, supplement p.19, equation (34-35);
# H and P should be computed before as in equations (28-29).
# Given p traits measured on n genotypes, Y needs to be a p x n matrix,
# and H and P need to be np x np

LL.REML.grad <- function(Y,P,Dk) {

  p <- nrow(Y)
  y <- matrix(Y)

  grad.G <- matrix(0,p,p)
  grad.E <- matrix(0,p,p)

  n <- ncol(Dk)

  for (j in 1:p) {
    for (i in 1:j) {
      ij.matrix        <- make.ij.matrix(i=i,j=j,p=p)
      grad.G[i,j]      <-  -0.5 * sum(diag(P %*% kronecker(Dk,ij.matrix)))      + 0.5 * t(y) %*% P %*% kronecker(Dk,     ij.matrix) %*% P %*% y
      grad.E[i,j]      <-  -0.5 * sum(diag(P %*% kronecker(diag(n),ij.matrix))) + 0.5 * t(y) %*% P %*% kronecker(diag(n),ij.matrix) %*% P %*% y
    }
  }

  diag(grad.G) <- diag(grad.G) / 2
  diag(grad.E) <- diag(grad.E) / 2

return(matrix(c(grad.G[upper.tri(grad.G,diag=T)],grad.E[upper.tri(grad.E,diag=T)])))
}


LL.hess <- function(Y,P,Dk,Hinv) {
  # Y=Yt
  p <- nrow(Y)

  n <- ncol(Dk)

  Q <- p * (p + 1) / 2

  y <- matrix(Y)

  hess.G <- matrix(0,Q,Q)

  hess.E <- matrix(0,Q,Q)

  hess.GE <- matrix(0,Q,Q)

  for (g2 in 1:Q) {
    for (g1 in 1:g2) {
      # g2 <- 1; g1 <- 1
      g1.indices <- which.i(i=g1,p=Q)
      g2.indices <- which.i(i=g2,p=Q)
      i          <- g1.indices[1]
      j          <- g1.indices[2]
      ip         <- g2.indices[1]
      jp         <- g2.indices[2]
      ij.matrix  <- make.ij.matrix(i=i,j=j,p=p,sym=T)
      ipjp.matrix<- make.ij.matrix(i=ip,j=jp,p=p,sym=T)

      G.trace            <- sum(diag(Hinv %*% kronecker(Dk,ij.matrix) %*% Hinv %*% kronecker(Dk,ipjp.matrix)))
      G.quad.form        <- t(y) %*% P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(Dk,ipjp.matrix) %*% P %*% y
      hess.G[g1,g2]      <-  (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * G.trace - G.quad.form)

      E.trace            <- sum(diag(Hinv %*% kronecker(diag(n),ij.matrix) %*% Hinv %*% kronecker(diag(n),ipjp.matrix)))
      E.quad.form        <- t(y) %*% P %*% kronecker(diag(n),ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix) %*% P %*% y
      hess.E[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * E.trace - E.quad.form)
      #hess.E[g1,g2]      <- E.quad.form#E.trace
    }
  }
  for (g2 in 1:Q) {
    for (g1 in 1:Q) {

      g1.indices <- which.i(i=g1,p=Q)
      g2.indices <- which.i(i=g2,p=Q)
      i          <- g1.indices[1]
      j          <- g1.indices[2]
      ip         <- g2.indices[1]
      jp         <- g2.indices[2]
      ij.matrix  <- make.ij.matrix(i=i,j=j,p=p)
      ipjp.matrix<- make.ij.matrix(i=ip,j=jp,p=p)

      GE.trace            <- sum(diag(Hinv %*% kronecker(Dk,ij.matrix) %*% Hinv %*% kronecker(diag(n),ipjp.matrix)))
      GE.quad.form        <- t(y) %*% P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix) %*% P %*% y
      hess.GE[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * GE.trace - GE.quad.form)
      #hess.GE[g1,g2]      <- GE.trace
    }
  }

  hess.G <- hess.G + t(hess.G); diag(hess.G) <- diag(hess.G) / 2
  hess.E <- hess.E + t(hess.E); diag(hess.E) <- diag(hess.E) / 2
  hessian <- rbind(cbind(hess.G,hess.GE),cbind(t(hess.GE),hess.E))
#return(list(hess.G,hess.E))
#return(c(hess.G[lower.tri(hess.G,diag=T)],hess.E[lower.tri(hess.E,diag=T)]))
return(hessian)
}

LL.REML.hess <- function(Y,P,Dk) {
  #Y=Yt
  p <- nrow(Y)

  n <- ncol(Dk)

  Q <- p * (p + 1) / 2

  y <- matrix(Y)

  hess.G <- matrix(0,Q,Q)

  hess.E <- matrix(0,Q,Q)

  hess.GE <- matrix(0,Q,Q)

  for (g2 in 1:Q) {
    for (g1 in 1:g2) {
      # g2 <- 1; g1 <- 1
      g1.indices <- which.i(i=g1,p=Q)
      g2.indices <- which.i(i=g2,p=Q)
      i          <- g1.indices[1]
      j          <- g1.indices[2]
      ip         <- g2.indices[1]
      jp         <- g2.indices[2]
      ij.matrix  <- make.ij.matrix(i=i,j=j,p=p)
      ipjp.matrix<- make.ij.matrix(i=ip,j=jp,p=p)

      G.trace            <- sum(diag(P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(Dk,ipjp.matrix)))
      G.quad.form        <- t(y) %*% P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(Dk,ipjp.matrix) %*% P %*% y
      hess.G[g1,g2]      <-  (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * G.trace - G.quad.form)

      E.trace            <- sum(diag(P %*% kronecker(diag(n),ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix)))
      E.quad.form        <- t(y) %*% P %*% kronecker(diag(n),ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix) %*% P %*% y
      hess.E[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * E.trace - E.quad.form)
    }
  }
  for (g2 in 1:Q) {
    for (g1 in 1:Q) {

      g1.indices <- which.i(i=g1,p=Q)
      g2.indices <- which.i(i=g2,p=Q)
      i          <- g1.indices[1]
      j          <- g1.indices[2]
      ip         <- g2.indices[1]
      jp         <- g2.indices[2]
      ij.matrix  <- make.ij.matrix(i=i,j=j,p=p)
      ipjp.matrix<- make.ij.matrix(i=ip,j=jp,p=p)

      GE.trace            <- sum(diag(P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix)))
      GE.quad.form        <- t(y) %*% P %*% kronecker(Dk,ij.matrix) %*% P %*% kronecker(diag(n),ipjp.matrix) %*% P %*% y
      hess.GE[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * GE.trace - GE.quad.form)
    }
  }

  hess.G <- hess.G + t(hess.G); diag(hess.G) <- diag(hess.G) / 2
  hess.E <- hess.E + t(hess.E); diag(hess.E) <- diag(hess.E) / 2
  hessian <- rbind(cbind(hess.G,hess.GE),cbind(t(hess.GE),hess.E))
#return(list(hess.G,hess.E))
#return(c(hess.G[lower.tri(hess.G,diag=T)],hess.E[lower.tri(hess.E,diag=T)]))
return(hessian)
}
