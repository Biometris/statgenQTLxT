#      Cm.new.output  <- updateC(s=Omega.2, rho=pen.Cm,V=V,sigma2.ind=sigma2.ind,PG=PG)

updateC <- function(Omega,rho,V,sigma2.ind.start,PG.start,tol.em=0.0001,
                    max.iter=200,print.progress=FALSE,penalize.diagonals=TRUE) {

#Omega=Omega.2;rho=pen.Cm;V=V;sigma2.ind.start=sigma2.ind;PG.start=PG; tol.em=0.0001;max.iter=300

  # cat('condi. number Omega')

  # When rank(Omega) = Q, W should be the Q x p matrix such that Omega = W^t W / Q
  W <- MatrixRoot(Omega)

  #W <- chol(Omega)

  Q <- nrow(W) # to do : - get the actual rank of W and truncate; now always Q=p
               #         - different notation: Q is confusing, given the Q.1 and Q.2 defined below

  # no fixed effects, so the sample covariance matrix Omega.2 = (1/Q) W^t W (not 1/(Q-...))
  W <- W * sqrt(Q)

  p <- ncol(W)

  B    <- NULL

  w <- eigen(diag(Q))

  Lambda.R <- diag(w$values)

  U <- w$vectors

  ############

  sigma2.ind <- sigma2.ind.start

  Cm  <- PG.start

  ############

  continue <- T

  iter <- 1

  E.log.lik <- 10000000000

  # to do : take the following outside the function
  Vinv <- ginv(V)

  VinvSqrt <- MatrixRoot(Vinv)

  VSqrt <- MatrixRoot(V)

  #############################################
  # EM

  while (continue & iter < max.iter) {

    # prevent that Cm become asymmetric because of numerical inaccuracies
    Cm <- (Cm + t(Cm)) / 2

    ########

    Dm.sqrt.inv <- sqrt(sigma2.ind) * VSqrt

    w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

    Q.1 <- w$vectors

    Lambda.1 <- w$values

    ##########

    Cm.sqrt.inv <- MatrixRoot(solve(Cm))

    w <- eigen(Cm.sqrt.inv %*% Vinv %*% Cm.sqrt.inv / sigma2.ind)

    Q.2 <- w$vectors

    Lambda.2 <- w$values

    ##########  In the preprint of Dahl et al (arxiv, version 6 dec. 2013),
    # part.1-part.4 are the quantities, on the bottom part of p.6, in this order
    # Each time we compute the right hand side of the equation
    # In dahl_etal_2013_debug.r we checked the left hand side(s) as well
    # Also S.1 and S.2 correspond to p. 6 of their preprint

    S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% W %*% VinvSqrt %*% Q.1 / sqrt(sigma2.ind))

    S.2     <- vec.inv.diag(x=Lambda.2,y=1/diag(Lambda.R))  *  (t(U) %*% W %*% MatrixRoot(Cm) %*% Q.2)

    part.1  <- Dm.sqrt.inv %*% Q.1 %*% diag(trace.p.diag.inv(x=Lambda.1,y=diag(Lambda.R))) %*% t(Q.1) %*% Dm.sqrt.inv

    part.2  <- Cm.sqrt.inv %*% Q.2 %*% diag(trace.p.diag.inv(x=Lambda.2,y=1/diag(Lambda.R))) %*% t(Q.2) %*% Cm.sqrt.inv

    part.3  <- (Cm.sqrt.inv %*% Q.2 %*% t(S.2)) %*% t(Cm.sqrt.inv %*% Q.2 %*% t(S.2))

    part.4  <- (Dm.sqrt.inv %*% Q.1 %*% t(S.1)) %*% Lambda.R %*% t(Dm.sqrt.inv %*% Q.1 %*% t(S.1))

    mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)


    ################
    # Compare with the 'naive' expressions:
    #Sigma <- ginv(kronecker((Vinv / sigma2.ind),diag(Q)) + kronecker(Cm,diag(Q))); M <- matrix(Sigma %*% (kronecker((Vinv / sigma2.ind),diag(Q))) %*% matrix(as.numeric(W)),ncol=p)
    #part1.check <- trace.p(Sigma,p,Q); part2.check <-  trace.p(kronecker(diag(p),diag(Q)) %*% Sigma,p,Q)
    #part3.check <- t(W - M) %*% (W - M); part4.check <- t(M) %*% diag(Q) %*% M
    #part1.check-part.1;part2.check-part.2;part3.check-part.3;part4.check-part.4
    #part.1;part.2;part.3;part.4
    #

    ###############

    Omega.1 <- (1 / Q) * part.3 + (1 / Q) * part.1

    Omega.2 <- (1 / Q) * part.2 + (1 / Q) * part.4

    ################   UPDATE :

    Cm.new <- glasso(s=Omega.2, rho=rho,penalize.diagonal=penalize.diagonals)$wi # ,penalize.diagonal=FALSE

    Cm.diff <- sum(abs(Cm.new - Cm))

    #########

    sigma2.ind.new <- (1 / p) * sum(diag(Vinv %*% Omega.1))

    sigma2.ind.diff <- abs(sigma2.ind.new - sigma2.ind)

    #########

    E.log.lik.old <- E.log.lik

    E.log.lik <- Q * log(det(Cm)) - Q * sum(diag(Cm %*% Omega.2)) + Q * log(det(Vinv / sigma2.ind)) - Q * sum(diag((Vinv / sigma2.ind) %*% Omega.1))

    if (print.progress) {
      cat('Iteration ',iter,' : ',Cm.diff,'  ',sigma2.ind.diff,'   ',E.log.lik,'\n')
    }

    Cm <- Cm.new

    sigma2.ind <- sigma2.ind.new

    iter <- iter + 1

    if (abs(E.log.lik - E.log.lik.old) < tol.em) {continue <- F}

  }

#cat('Inner loop: ',sum(abs(part.1)),' ',sum(abs(part.2)),' ',sum(abs(part.3)),' ',sum(abs(part.4)),'\n')

return(list(PG=Cm,sigma2.ind=sigma2.ind,n.iter=iter))
}