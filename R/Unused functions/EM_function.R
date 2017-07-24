EM_function <- function(Y,K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,
                        tol.em=0.0001,max.iter.em=300,
                        pen.Cm=0,pen.Dm=0,
                        Cm.start=NULL,Dm.start=NULL,
                        L2.penalty=F,
                        L2.target.Cm=0.5 * diag(ncol(Y)),
                        L2.target.Dm=0.5 * diag(ncol(Y)),
                        stop.if.not.monotone=FALSE) {
# implements the EM-algorithm of Dahl et al 2013, arxiv

# X default used to be : X=data.frame()

# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed
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

#Y=Y;K=K;X=matrix(rep(1,nrow(K)));Dm.is.diagonal=F;tol.em=0.0001;max.iter.em=300;pen.Cm=0;pen.Dm=0;Cm.start=NULL;Dm.start=NULL;L2.penalty=F;

# Y=Y[nonmissing.accessions,2:3];K=K[nonmissing.accessions,nonmissing.accessions];X=matrix(rep(1,length(nonmissing.accessions)));Dm.is.diagonal=T;pen.Cm=0;pen.Dm=0;max.iter.em=2000;tol.em=0.00001;Cm.start=NULL;Dm.start=NULL;L2.penalty=F
  require(MASS)


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
    Cm <- diag(apply(Y,2,var)) / 2#diag(p) / 2#ginv((cor(Y) + diag(p)) / 4) #0.5 * cor(Y)
  } else {
    Cm <- Cm.start
  }

  if (is.null(Dm.start)) {
    Dm <- diag(apply(Y,2,var)) / 2 #diag(p) / 2
  } else {
    Dm <- Dm.start
  }

  continue <- T

  iter <- 1

  E.log.lik <-  - 10000000000

  #############################################
  # EM

  while (continue & iter < max.iter.em) {

    # prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- (Cm + t(Cm)) / 2

    Dm <- (Dm + t(Dm)) / 2

    ########

    #Dm.sqrt.inv <- sqrtm(ginv(Dm))
    Dm.sqrt.inv <- MatrixRoot(ginv(Dm))
    # before 15-6-2015
    #Dm.sqrt.inv <- MatrixRoot(solve(Dm))

    w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

    Q.1 <- w$vectors

    Lambda.1 <- w$values

    ##########

    #Cm.sqrt.inv <- sqrtm(ginv(Cm))
    Cm.sqrt.inv <- MatrixRoot(ginv(Cm))
    # before 15-6-2015
    #Cm.sqrt.inv <- MatrixRoot(solve(Cm))


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

    ################
    # Compare with the 'naive' expressions:
    #Sigma <- ginv(kronecker(Dm,diag(n)) + kronecker(Cm,R)); M <- matrix(Sigma %*% (kronecker(Dm,diag(n))) %*% matrix(as.numeric(Y)),ncol=p)
    #part1.check <- trace.p(Sigma,p,n); part2.check <-  trace.p(kronecker(diag(p),R) %*% Sigma,p,n)
    #part3.check <- t(Y - M) %*% (Y - M); part4.check <- t(M) %*% R %*% M
    #part1.check;part2.check;part3.check;part4.check
    #part.1;part.2;part.3;part.4
    #
    ###############

    mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)

    #mu <- matrix(rep(0,n*p),ncol=p)

    if (nc > 0) {
      B  <- XtXinvXt %*% (Y - mu)
    }

    Omega.1 <- (1 / n) * part.3 + (1 / n) * part.1

    Omega.2 <- (1 / n) * part.2 + (1 / n) * part.4


    ################

    if (L2.penalty) {
      Cm.new <- ridgeS(S=Omega.2, lambda=pen.Cm, target = L2.target.Cm)
    } else {
      Cm.new <- glasso(s=Omega.2, rho=pen.Cm)$wi
    }

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

    E.log.lik <- log(det(Cm)) - sum(diag(Cm %*% Omega.2)) + log(det(Dm)) - sum(diag(Dm %*% Omega.1))

    pen <- pen.Cm * sum(abs(Cm)) + pen.Dm * sum(abs(Dm))

    cat('Iteration ',iter,' : ',Cm.diff,'  ',Dm.diff,'   ',E.log.lik - pen,'\n')

    #LL <- 0
    #for (j in 1:p) {
    #  Vj  <- K / diag(Cm)[j] + diag(n) / diag(Dm)[j]
    #  yj  <- matrix(Y[,j])
    #  LLj <- -0.5 * log(det(Vj)) -0.5 * as.numeric(t(yj) %*% solve(Vj) %*% yj)
    #  LL  <- LL + LLj
    #}
    #cat(LL,'\n')

    Cm <- Cm.new

    Dm <- Dm.new

    iter <- iter + 1

    #if (Cm.diff < tol & Dm.diff < tol) {continue <- F}
    if (abs(E.log.lik - E.log.lik.old) < tol.em) {continue <- F}

    if (stop.if.not.monotone) {if (E.log.lik < E.log.lik.old) {stop('Error: decreasing LL')}}

    # (Y - X %*% B)

  }

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype',1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait',1:p)}

  pred.frame <- data.frame(trait=rep(colnames(Y),each=n),genotype=rep(rownames(Y),p),predicted=as.numeric(mu))

return(list(Cm=Cm,Dm=Dm,B=B,pred=pred.frame,log.lik=E.log.lik,tol.em=tol.em,converged=(!continue),
            pen.Cm=pen.Cm,pen.Dm=pen.Dm))
}
#test <- EM_function(Y=Y,K=K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,pen.Cm=2,pen.Dm=2)

############################################
