EM_function_with_indices <- function(Y,K,V,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,
                        penalize.diagonals=TRUE,
                        tol.em=0.0001,max.iter.em=300,
                        pen.Cm=0,pen.Dm=0,
                        Cm.start=NULL,Dm.start=NULL,
                        L2.penalty=F,
                        L2.target.Cm=0.5 * diag(ncol(Y)),
                        L2.target.Dm=0.5 * diag(ncol(Y))) {
# implements the EM-algorithm of Dahl et al 2013, arxiv; with an additional environmental matrix describing covariance described by indices

# Y=Ysim;X=data.frame();Dm.is.diagonal=T; pen.Cm=0;pen.Dm=0;L2.penalty=F; Cm.start=NULL;Dm.start=NULL

# X default used to be : X=data.frame()

# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# V            : the 'environmental' kinship matrix
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed
# penalize.diagonals : if FALSE, do not penalize the diagonals of Cm and Dm, when using glasso
# Dm.is.diagonal:
# tol.em       : tolerance in the EM-algorithm : stop when the increase in log-lik. is less than tol
# max.iter.em  : maximum number of iterations in the EM algorithm
# pen.Cm       : penalty on the precision matrix Cm of genetic correlations
# pen.Dm       : penalty on the precision matrix Dm of environmental correlations
# Cm.start     : starting values for Cm, as p x p matrix
# Dm.start     : starting values for Dm, as p x p matrix
# L2.penalty   : if TRUE, ridge estimation is performed (van Wieringen and Peeters, 2014)
#                shrinking to L2.target.Cm(Dm)
#               If FALSE, an L1 penalty is used as in Dahl et al (2013)

  if (L2.penalty) {require(rags2ridges)} else {require(glasso)}

  # also check : missing data

  stopifnot(nrow(Y)==nrow(K))
  stopifnot(ncol(K)==nrow(K))

  nc    <- ncol(X)
  if (nc > 0) {stopifnot(nrow(X)==nrow(K))}

  n <- ncol(K)

  p <- ncol(Y)

  if (nc > 0) {
    B     <- matrix(0,nc,p) # the c x p matrix of coefficients (p traits)
    XtXinvXt <- solve(t(X) %*% X) %*% t(X)
  } else {
    B    <- NULL
  }

  R <- solve(K)

  w <- eigen(R)

  Lambda.R <- diag(w$values)

  U <- w$vectors

  if (is.null(Cm.start)) {
    Cm <- ginv((cor(Y) + diag(p)) / 4) #0.5 * cor(Y)
  } else {
    Cm <- Cm.start
  }

  ############
  # the model is Cm^{-1} = PG^{-1} + sigma_ind * V =
  # Given a starting value for Cm, set starting values for sigma2.ind and PG

  InvCm <- ginv(Cm)

  sigma2.ind <- sum(diag(InvCm)) / sum(diag(V))

  while (!is.positive.definite(InvCm - sigma2.ind * V)) {sigma2.ind <- 0.5 * sigma2.ind}

  PG  <- ginv(InvCm - sigma2.ind * V)

  ############


  if (is.null(Dm.start)) {
    Dm <- diag(p) / 2
  } else {
    Dm <- Dm.start
  }

  continue <- T

  iter <- 1

  E.log.lik <- 10000000000

  #############################################
  # EM

  while (continue & iter < max.iter.em) {

    # prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- (Cm + t(Cm)) / 2

    Dm <- (Dm + t(Dm)) / 2

    ########

    Dm.sqrt.inv <- MatrixRoot(solve(Dm))

    w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

    Q.1 <- w$vectors

    Lambda.1 <- w$values

    ##########

    Cm.sqrt.inv <- MatrixRoot(solve(Cm))

    ################################### Some lines for debugging ##########
    #print(sum(abs(Cm-t(Cm))));cat('\n')
    #print(sum(abs(Dm-t(Dm))));cat('\n')
    #if (iter > 198) {
    #  print(Cm.sqrt.inv)
    #  cat('\n')
    #  print(Dm.sqrt.inv)
    #  cat('\n')
    #  #print(sum(abs(Cm-t(Cm))))#eigen(solve(Cm))$values)#eigen(Cm)$values)
    #  cat('\n')
    #  print(Dm)
    #  cat('\n')
    #}
    #######################################################################

    w <- eigen(Cm.sqrt.inv %*% Dm %*% Cm.sqrt.inv)

    Q.2 <- w$vectors

    Lambda.2 <- w$values

    ##########  In the preprint of Dahl et al (arxiv, version 6 dec. 2013),
    # part.1-part.4 are the quantities, on the bottom part of p.6, in this order
    # Each time we compute the right hand side of the equation
    # In dahl_etal_2013_debug.r we checked the left hand side(s) as well
    # Also S.1 and S.2 correspond to p. 6 of their preprint

    if (nc > 0) {

      S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% (Y - X %*% B) %*% MatrixRoot(Dm) %*% Q.1)

      S.2     <- vec.inv.diag(x=Lambda.2,y=1/diag(Lambda.R))  *  (t(U) %*% (Y - X %*% B) %*% MatrixRoot(Cm) %*% Q.2)

    } else {

      S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% Y %*% MatrixRoot(Dm) %*% Q.1)

      S.2     <- vec.inv.diag(x=Lambda.2,y=1/diag(Lambda.R))  *  (t(U) %*% Y %*% MatrixRoot(Cm) %*% Q.2)

    }

    part.1  <- Dm.sqrt.inv %*% Q.1 %*% diag(trace.p.diag.inv(x=Lambda.1,y=diag(Lambda.R))) %*% t(Q.1) %*% Dm.sqrt.inv

    part.2  <- Cm.sqrt.inv %*% Q.2 %*% diag(trace.p.diag.inv(x=Lambda.2,y=1/diag(Lambda.R))) %*% t(Q.2) %*% Cm.sqrt.inv

    part.3  <- (Cm.sqrt.inv %*% Q.2 %*% t(S.2)) %*% t(Cm.sqrt.inv %*% Q.2 %*% t(S.2))

    part.4  <- (Dm.sqrt.inv %*% Q.1 %*% t(S.1)) %*% Lambda.R %*% t(Dm.sqrt.inv %*% Q.1 %*% t(S.1))

    ###############

    mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)

    if (nc > 0) {
      B  <- XtXinvXt %*% (Y - mu)
    }


    ################
    # Compare with the 'naive' expressions:
    #Sigma <- ginv(kronecker(Dm,diag(n)) + kronecker(Cm,R)); M <- matrix(Sigma %*% (kronecker(Dm,diag(n))) %*% matrix(as.numeric(Y)),ncol=p)
    #part1.check <- trace.p(Sigma,p,n); part2.check <-  trace.p(kronecker(diag(p),R) %*% Sigma,p,n)
    #part3.check <- t(Y - M) %*% (Y - M); part4.check <- t(M) %*% R %*% M
    #part1.check-part.1;part2.check-part.2;part3.check-part.3;part4.check-part.4
    #part.1;part.2;part.3;part.4
    ##############################

    Omega.1 <- (1 / n) * part.3 + (1 / n) * part.1

    Omega.2 <- (1 / n) * part.2 + (1 / n) * part.4


    ################

    if (L2.penalty) {
      stop('The L2 penalty has not yet been implemented in this function')
      Cm.new <- ridgeS(S=Omega.2, lambda=pen.Cm, target = L2.target.Cm)
    } else {
      #Cm.new        <- glasso(s=Omega.2, rho=pen.Cm)$wi
      Cm.new.output  <- updateC(Omega=Omega.2, rho=pen.Cm,V=V,
                                sigma2.ind.start=sigma2.ind,PG.start=PG,
                                penalize.diagonals=penalize.diagonals)

      sigma2.ind.new <- Cm.new.output$sigma2.ind
      PG.new         <- Cm.new.output$PG
      Cm.new         <- ginv(ginv(PG.new) + sigma2.ind.new * V)

    }

    ##################

    if (Dm.is.diagonal) {
      tau    <- p / sum(diag(Omega.1))
      Dm.new <- tau * diag(p)
    } else {
      if (L2.penalty) {
        Dm.new <- ridgeS(S=Omega.1, lambda=pen.Dm, target = L2.target.Dm)
      } else {
        Dm.new <- glasso(s=Omega.1, rho=pen.Dm)$wi
      }
    }

    Cm.diff <- sum(abs(Cm.new - Cm))

    Dm.diff <- sum(abs(Dm.new - Dm))

    E.log.lik.old <- E.log.lik

    E.log.lik <- n * log(det(Cm)) - n * sum(diag(Cm %*% Omega.2)) + n * log(det(Dm)) - n * sum(diag(Dm %*% Omega.1))

    cat('Iteration ',iter,' : ',Cm.diff,'  ',Dm.diff,'   ',sigma2.ind.new,'    ',Cm.new.output$n.iter,'   ',E.log.lik,'\n')

    cat('smallest eigenvalue Omega.1: ',max(eigen(Omega.1)$values)/min(eigen(Omega.1)$values),'\n')
    cat('smallest eigenvalue Omega.2: ',max(eigen(Omega.2)$values)/min(eigen(Omega.2)$values),'\n')
    cat('smallest eigenvalue Cm: ',min(eigen(Cm)$values),'\n')
    cat('smallest eigenvalue Dm: ',min(eigen(Dm)$values),'\n')
    cat('L1 norm Cm:           : ',sum(abs(Cm)),'\n')
    cat('L1 norm Dm:           : ',sum(abs(Dm)),'\n')
    cat('Outer loop: ',sum(abs(part.1)),' ',sum(abs(part.2)),' ',sum(abs(part.3)),' ',sum(abs(part.4)),'\n')

    #############

    Cm <- Cm.new

    Dm <- Dm.new

    sigma2.ind <- sigma2.ind.new

    PG <- PG.new

    #############

    iter <- iter + 1

    #if (Cm.diff < tol & Dm.diff < tol) {continue <- F}
    if (abs(E.log.lik - E.log.lik.old) < tol.em) {continue <- F}

  }

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype',1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait',1:p)}

  pred.frame <- data.frame(trait=rep(colnames(Y),each=n),genotype=rep(rownames(Y),p),predicted=as.numeric(mu))

return(list(Cm=Cm,sigma2.ind=sigma2.ind,PG=PG,Dm=Dm,B=B,pred=pred.frame,log.lik=E.log.lik,tol.em=tol.em,converged=(!continue),
            pen.Cm=pen.Cm,pen.Dm=pen.Dm))
}
#test <- EM_function(Y=Y,K=K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,pen.Cm=2,pen.Dm=2)

