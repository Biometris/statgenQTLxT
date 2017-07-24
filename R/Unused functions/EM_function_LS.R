# LS = low rank + sparse
#
EM_function_LS <- function(Y,K,X=matrix(rep(1,nrow(K))),
                        Dm.is.diagonal=F,
                        penalize.diagonal.G=TRUE,
                        penalize.diagonal.E=TRUE,
                        tol.em=0.0001,
                        max.iter.em=300,
                        pen.Cm=0,
                        pen.Dm=0,
                        Cm.start=NULL,
                        Dm.start=NULL,
                        m.G=1,
                        m.E=1,
                        prediction=FALSE) {

#Y=Ysim;K=K;X=data.frame();Dm.is.diagonal=F;pen.Cm=100;pen.Dm=1000;max.iter.em=max.iter.em;tol.em=tol.em;penalize.diagonal.G=FALSE;penalize.diagonal.E=FALSE;Cm.start=NULL;Dm.start=NULL;m.G=2;m.E=1

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

  require(glasso)

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

  ############
  # check if m.G and m.E have sensible values

  if (m.G!=round(m.G)) {stop("m.G needs to be integer")}
  if (m.E!=round(m.E)) {stop("m.E needs to be integer")}

  if (m.G < 0) {stop("m.G cannot be negative")}
  if (m.E < 0) {stop("m.E cannot be negative")}

  if (m.G >= p) {stop("m.G needs to be smaller than the number of traits or environments")}
  if (m.E >= p) {stop("m.E needs to be smaller than the number of traits or environments")}

  #############

  if (m.E > 0 & Dm.is.diagonal==TRUE) {
    Dm.is.diagonal <- FALSE
    cat('Warning: m.E is larger than zero, hence Dm.is.diagonal is set to FALSE \n')
  }

  ###############

  if (is.null(Cm.start)) {
    Cm <- ginv((cor(Y) + diag(p)) / 4)
  } else {
    Cm <- Cm.start
  }


  if (is.null(Dm.start)) {
    if (Dm.is.diagonal) {
      Dm <- diag(p) / 2
    } else {
      Dm <- ginv((cor(Y) + diag(p)) / 4)
    }
  } else {
    Dm <- Dm.start
  }


  ############
  # the model is Cm^{-1} = P^{-1} + W W^t =
  # Given a starting value for Cm, set starting values for P and W

  if (m.G > 0) {

    InvCm      <- ginv(Cm)
    eigen.temp <- eigen(InvCm)
    U.G        <- as.matrix(eigen.temp$vectors[,1:m.G])
    psi.G      <- mean(eigen.temp$values[-(1:m.G)])
    if (m.G > 1) {
      root.Lambda.G   <- MatrixRoot(diag(eigen.temp$values[1:m.G] - rep(psi.G,m.G)))
    }
    if (m.G == 1) {
      root.Lambda.G   <- matrix(sqrt(eigen.temp$values[1:m.G] - rep(psi.G,m.G)))
    }
    W.G        <- U.G %*% root.Lambda.G
    if (is.positive.definite(InvCm - W.G %*% t(W.G))) {
      P.G        <- solve(InvCm - W.G %*% t(W.G))
    } else {
      P.G <- diag(p) / psi.G
    }

    # check
    # eigen(InvCm - W.G %*% t(W.G)) ; W.G %*% t(W.G)
  } else {

    W.G <- NULL
    P.G <- NULL

  }

  if (m.E > 0) {

    InvDm      <- ginv(Dm)
    eigen.temp <- eigen(InvDm)
    U.E        <- as.matrix(eigen.temp$vectors[,1:m.E])
    psi.E      <- mean(eigen.temp$values[-(1:m.E)])
    if (m.E > 1) {
      root.Lambda.E   <- MatrixRoot(diag(eigen.temp$values[1:m.E] - rep(psi.E,m.E)))
    }
    if (m.E == 1) {
      root.Lambda.E   <- matrix(sqrt(eigen.temp$values[1:m.E] - rep(psi.E,m.E)))
    }
    W.E        <- U.E %*% root.Lambda.E
    if (is.positive.definite(InvCm - W.E %*% t(W.E))) {
      P.E        <- solve(InvCm - W.E %*% t(W.E))
    } else {
      P.E <- diag(p) / psi.E
    }
    # check
    # eigen(InvCm - W.E %*% t(W.E)) ; W.E %*% t(W.E)
  } else {

    W.E <- NULL
    P.E <- NULL

  }

  ############

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

    Dm.sqrt.inv <- MatrixRoot(ginv(Dm)) #MatrixRoot(solve(Dm))

    w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

    Q.1 <- w$vectors

    Lambda.1 <- w$values

    ##########

    Cm.sqrt.inv <- MatrixRoot(ginv(Cm)) # MatrixRoot(solve(Cm))

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

    if (prediction) {
      mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)
    } else {
      mu <- matrix(rep(0,n*p),ncol=p)
    }

    #mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)

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

    ################  update C

    if (m.G == 0) {

      Cm.new <- glasso(s=Omega.2, rho=pen.Cm)$wi

    } else {

      # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- MatrixRoot(Omega.2)
      A <- A * sqrt(nrow(A))

      Cm.new.output  <- update_LS(Y=A, rho=pen.Cm,
                                  W.start=W.G,
                                  P.start=P.G,
                                  penalize.diagonal=penalize.diagonal.G)

      W.G.new         <- Cm.new.output$W
      P.G.new         <- Cm.new.output$P

      Cm.new          <- ginv(ginv(P.G.new) + W.G.new %*% t(W.G.new))

    }


    ################## update D

    if (Dm.is.diagonal) {

      tau    <- p / sum(diag(Omega.1))
      Dm.new <- tau * diag(p)

    } else {

      if (m.E == 0) {

        Dm.new <- glasso(s=Omega.1, rho=pen.Dm)$wi

      } else {

        # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
        A <- MatrixRoot(Omega.1)
        A <- A * sqrt(nrow(A))

        Dm.new.output  <- update_LS(Y=A, rho=pen.Dm,
                                    W.start=W.E,
                                    P.start=P.E,
                                    penalize.diagonal=penalize.diagonal.E)

        W.E.new         <- Dm.new.output$W
        P.E.new         <- Dm.new.output$P

        Dm.new          <- ginv(ginv(P.E.new) + W.E.new %*% t(W.E.new))
      }
    }

    ####################

    Cm.diff <- sum(abs(Cm.new - Cm))

    Dm.diff <- sum(abs(Dm.new - Dm))

    E.log.lik.old <- E.log.lik

    E.log.lik <- n * log(det(Cm)) - n * sum(diag(Cm %*% Omega.2)) + n * log(det(Dm)) - n * sum(diag(Dm %*% Omega.1))

    cat('Iteration ',iter,' : ',Cm.diff,'  ',Dm.diff,'    ',E.log.lik,'\n') # Cm.new.output$n.iter,'   '

    #############

    Cm <- Cm.new

    Dm <- Dm.new

    #############

    iter <- iter + 1

    #if (Cm.diff < tol & Dm.diff < tol) {continue <- F}
    if (abs(E.log.lik - E.log.lik.old) < tol.em) {continue <- F}

  }

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype',1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait',1:p)}

  pred.frame <- data.frame(trait=rep(colnames(Y),each=n),genotype=rep(rownames(Y),p),predicted=as.numeric(mu))

return(list(Cm=Cm,Dm=Dm,B=B,
            W.G=W.G,W.E=W.E,
            P.G=P.G,P.E=P.E,
            pred=pred.frame,log.lik=E.log.lik,
            tol.em=tol.em,converged=(!continue),
            pen.Cm=pen.Cm,pen.Dm=pen.Dm))
}
#test <- EM_function(Y=Y,K=K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,pen.Cm=2,pen.Dm=2)

