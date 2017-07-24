# FA = factor analytic ; a variation on the more general function EM_function_LS (low rank + sparse)
#
EM_function_FA <- function(Y,K,X=matrix(rep(1,nrow(K))),
                        Cm.het=F,Dm.het=F,
                        tol.em=0.0001,
                        max.iter.em=300,
                        Cm.start=NULL,
                        Dm.start=NULL,
                        m.G=1,
                        m.E=1,
                        max.diag=1e04,
                        prediction=TRUE,
                        stopifdecreasing=FALSE,
                        compute.log.lik=FALSE) {


# Cm.start=NULL; Dm.start=NULL; m.G=1; m.E=1; Cm.het=T; Dm.het=T;compute.log.lik=TRUE; stopifdecreasing=FALSE

# Y=Y;K=K;X=X;max.iter.em=max.iter.em;tol.em=tol.em;Cm.start=NULL;Dm.start=NULL;m.G=m.G;m.E=m.E;Cm.het=TRUE;Dm.het=TRUE;max.diag=100; prediction=TRUE;stopifdecreasing=FALSE;compute.log.lik=T

# X default used to be : X=data.frame()

# Y            : the n x p matrix of phenotypic observations (n individuals, p traits)
#                without missing values; NOT transformed
# K            : the n x n kinship matrix; NOT transformed
# X            : the n x c design matrix  (n individuals, c covariates), can be data.frame()
#                NOT transformed
# Dm.is.diagonal:
# tol.em       : tolerance in the EM-algorithm : stop when the increase in log-lik. is less than tol
# max.iter.em  : maximum number of iterations in the EM algorithm
# Cm.start     : starting values for Cm, as p x p matrix
# Dm.start     : starting values for Dm, as p x p matrix

  # also check : missing data

  Y <- as.matrix(Y)

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

  #R <- ginv(K)
  R <- Ginv(K)


  w <- eigen(R)

  #eigen.values.K <- (1/eigen(R)$values,dec=T)

  w.K <- eigen(K)
  Dk  <- diag(w$values)
  Uk  <- w$vectors

  Lambda.R <- as.matrix(diag(w$values))

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

  #if (m.E > 0 & Dm.is.diagonal==TRUE) {
  #  Dm.is.diagonal <- FALSE
  #  cat('Warning: m.E is larger than zero, hence Dm.is.diagonal is set to FALSE \n')
  #}

  ###############

  if (is.null(Cm.start)) {
    if (m.G==0) {
      Cm <- 2 * as.matrix(diag(p))
    } else {
      #Cm <- ginv((cor(Y) + diag(p)) / 4)
      Cm <- Ginv((cor(Y) + diag(p)) / 4)
    }
  } else {
    Cm <- Cm.start
  }

  if (is.null(Dm.start)) {
    if (m.E==0) {
      Dm <- 2 * as.matrix(diag(p))
    } else {
      #Dm <- ginv((cor(Y) + diag(p)) / 4)
      Dm <- Ginv((cor(Y) + diag(p)) / 4)
    }
  } else {
    Dm <- Dm.start
  }

  ############
  # the model is Cm^{-1} = P^{-1} + W W^t =
  # Given a starting value for Cm, set starting values for P and W

  if (m.G > 0) {

    #InvCm      <- ginv(Cm)
    InvCm      <- Ginv(Cm)
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
    P.G        <- diag(p) / psi.G

  } else {

    W.G <- NULL
    P.G <- NULL

  }

  if (m.E > 0) {

    #InvDm      <- ginv(Dm)
    InvDm      <- Ginv(Dm)
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
    P.E        <- diag(p) / psi.E

  } else {

    W.E <- NULL
    P.E <- NULL

  }

  ############

  continue <- TRUE

  decreased <- FALSE

  iter <- 1

  E.log.lik.Cm <- -10000000000

  E.log.lik.Dm <- -10000000000

  E.log.lik    <- E.log.lik.Cm + E.log.lik.Dm

  mu <- matrix(rep(0,n*p),ncol=p)

  #############################################
  # EM

  while (continue & iter < max.iter.em) {

    # prevent that Cm, Dm become asymmetric because of numerical inaccuracies
    Cm <- as.matrix((Cm + t(Cm))) / 2

    Dm <- as.matrix((Dm + t(Dm))) / 2

    ########

    Dm.sqrt.inv <- MatrixRoot(Ginv(Dm))
    #Dm.sqrt.inv <- try(MatrixRoot(ginv(Dm)), silent=TRUE)
    #if (class(Dm.sqrt.inv)=="try-error") {
    #  Dm.sqrt.inv <- MatrixRoot(solve(Dm))
    #}

    w <- eigen(Dm.sqrt.inv %*% Cm %*% Dm.sqrt.inv)

    Q.1 <- w$vectors

    Lambda.1 <- w$values

    ##########

    Cm.sqrt.inv <- MatrixRoot(Ginv(Cm)) # MatrixRoot(solve(Cm))


    #Cm.sqrt.inv <- try(MatrixRoot(ginv(Cm)), silent=TRUE)
    #if (class(Cm.sqrt.inv)=="try-error") {
    #  Cm.sqrt.inv <- MatrixRoot(solve(Cm))
    #}

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

      tUYminXb <- t(U) %*% (Y - X %*% B)

      S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (tUYminXb %*% (MatrixRoot(Dm) %*% Q.1))

      S.2     <- vec.inv.diag(x=Lambda.2,y=1/diag(Lambda.R))  *  (tUYminXb %*% (MatrixRoot(Cm) %*% Q.2))

    } else {

      S.1     <- vec.inv.diag(x=Lambda.1,y=diag(Lambda.R))    *  (t(U) %*% Y %*% MatrixRoot(Dm) %*% Q.1)

      S.2     <- vec.inv.diag(x=Lambda.2,y=1/diag(Lambda.R))  *  (t(U) %*% Y %*% MatrixRoot(Cm) %*% Q.2)

    }


    if (p > 1) {

      part.1  <- Dm.sqrt.inv %*% Q.1 %*% (diag(trace.p.diag.inv(x=Lambda.1,y=diag(Lambda.R))) %*% t(Q.1) %*% Dm.sqrt.inv)

      part.2  <- Cm.sqrt.inv %*% Q.2 %*% (diag(trace.p.diag.inv(x=Lambda.2,y=1/diag(Lambda.R))) %*% t(Q.2) %*% Cm.sqrt.inv)

    } else {

      part.1  <- Dm.sqrt.inv %*% Q.1 %*% (matrix(trace.p.diag.inv(x=Lambda.1,y=diag(Lambda.R))) %*% (t(Q.1) %*% Dm.sqrt.inv))

      part.2  <- Cm.sqrt.inv %*% Q.2 %*% (matrix(trace.p.diag.inv(x=Lambda.2,y=1/diag(Lambda.R))) %*% (t(Q.2) %*% Cm.sqrt.inv))

    }

    part.3  <- Cm.sqrt.inv %*% Q.2 %*% t(S.2) %*% t(Cm.sqrt.inv %*% (Q.2 %*% t(S.2)))

    part.4  <- Dm.sqrt.inv %*% Q.1 %*% t(S.1) %*% Lambda.R %*% t(Dm.sqrt.inv %*% (Q.1 %*% t(S.1)))

    ###############

#    if (prediction) {
#      #mu <- matrix(kronecker(Dm.sqrt.inv %*% Q.1, U) %*% matrix(S.1),ncol=p)
#      mu <- matrix(U %*% S.1 %*% t(Dm.sqrt.inv %*% Q.1),ncol=p)
#    } else {
#      mu <- matrix(rep(0,n*p),ncol=p)
#    }

    if (nc > 0) {
      mu <- matrix(U %*% S.1 %*% t(Dm.sqrt.inv %*% Q.1),ncol=p)
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

    Omega.1 <- as.matrix((1 / n) * part.3 + (1 / n) * part.1)

    Omega.2 <- as.matrix((1 / n) * part.2 + (1 / n) * part.4)

    ################  update C

    if (m.G == 0) {

      # Recall that the model is Cm^{-1} = P^{-1} + W W^t.
      #
      # when m.G==0, W=0 and Cm=P
      if (p > 1) {
        if (Cm.het) {
          P.G.new <- diag(pmin(max.diag,1 / diag(Omega.2)))
        } else {
          tau    <- min(max.diag,p / sum(diag(Omega.2)))
          P.G.new <- tau * diag(p)
        }
      } else {
        P.G.new <- matrix(1 / as.numeric(Omega.2))
        #P.G.new <- matrix(pmin(max.diag,1 / as.numeric(Omega.2)))
      }


      W.G.new <- NULL
      Cm.new  <- P.G.new

    } else {

      # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- MatrixRoot(Omega.2)

      A <- A * sqrt(nrow(A))

      if (Cm.het==FALSE) {

        #Cm.new.output   <- update_FA_homogeneous_var(Y=A,m=m.G)
        Cm.new.output   <- update_FA_homogeneous_var(S=Omega.2,m=m.G)
        Cm.new.output$P <- diag(p) / Cm.new.output$sigma2

      } else {

        Cm.new.output  <- update_FA(Y=A,
                                    W.start=W.G,
                                    P.start=P.G,
                                    het.var=Cm.het,
                                    max.diag=max.diag)
      }

      W.G.new         <- Cm.new.output$W
      P.G.new         <- Cm.new.output$P

      #Cm.new          <- ginv(ginv(P.G.new) + W.G.new %*% t(W.G.new))
      Cm.new          <- Ginv(Ginv(P.G.new) + W.G.new %*% t(W.G.new))
    }

    ################  update D

    if (m.E == 0) {


      # Recall that the model is Dm^{-1} = P^{-1} + W W^t.
      #
      # when m.E==0, W=0 and Cm=P

      if (p > 1) {

        if (Dm.het) {
          P.E.new <- diag(pmin(max.diag,1 / diag(Omega.1)))
          #P.E.new <- diag(1 / diag(Omega.1))
        } else {
          tau    <- min(max.diag,p / sum(diag(Omega.1)))
          P.E.new <- tau * diag(p)
        }

      } else {

        P.E.new <- matrix(1 / as.numeric(Omega.1))
        #P.E.new <- matrix(pmin(max.diag,1 / as.numeric(Omega.1)))
      }



      W.E.new <- NULL
      Dm.new  <- P.E.new

    } else {

      # When rank(Omega) = Q, A should be the Q x p matrix such that Omega = A^t A / Q
      A <- MatrixRoot(Omega.1)
      A <- A * sqrt(nrow(A))

      if (Dm.het==FALSE) {

        #Dm.new.output  <- update_FA_homogeneous_var(Y=A,m=m.E)
        Dm.new.output  <- update_FA_homogeneous_var(S=Omega.1,m=m.E)
        Dm.new.output$P <- diag(p) / Dm.new.output$sigma2

      } else {

        Dm.new.output  <- update_FA(Y=A,
                                    W.start=W.E,
                                    P.start=P.E,
                                    het.var=Dm.het,
                                    max.diag=max.diag)

        #Dm.new.output$P <- diag(p) / Dm.new.output$sigma2

      }
      W.E.new         <- Dm.new.output$W
      P.E.new         <- Dm.new.output$P

      #Dm.new          <- ginv(ginv(P.E.new) + W.E.new %*% t(W.E.new))
      Dm.new          <- Ginv(Ginv(P.E.new) + W.E.new %*% t(W.E.new))
    }

    ####################

    Cm.diff <- sum(abs(Cm.new - Cm))

    Dm.diff <- sum(abs(Dm.new - Dm))

    #######################

    log.lik <- NA

    # X is the c x n covariate matrix, c being the number of covariates and
    #      n being the number of genotypes
    #      c has to be at least one (typically an intercept)
    # Y is the p x n matrix of observed phenotypes, on p traits or environments
    #
    # No missing values are allowed in any of X,Y and x
    #
    # It is assumed that X,Y have already been rotated by Uk, where Uk is such that
    # the kinship matrix equals K = Uk %*% Dk %*% t(Uk)
    # (R: w <- eigen(K); Dk <- diag(w$values); Uk <- w$vectors)
    # Next, the original X, Y are post-multiplied by Uk, e.g. Y <- Y %*% Uk;

    #######################

    E.log.lik.old.Cm <- E.log.lik.Cm

    E.log.lik.old.Dm <- E.log.lik.Dm

    E.log.lik.old    <- E.log.lik.old.Cm + E.log.lik.old.Dm

    #E.log.lik <- n * log(det(Cm)) - n * sum(diag(Cm %*% Omega.2)) + n * log(det(Dm)) - n * sum(diag(Dm %*% Omega.1))

    E.log.lik.Cm     <- n * log(det(Cm)) - n * sum(diag(Cm %*% Omega.2))

    E.log.lik.Dm     <- n * log(det(Dm)) - n * sum(diag(Dm %*% Omega.1))

    #if (E.log.lik.Cm < E.log.lik.old.Cm) {stop('Error in updating Cm')}
    #if (E.log.lik.Dm < E.log.lik.old.Dm) {stop('Error in updating Dm')}

    E.log.lik        <- E.log.lik.Cm + E.log.lik.Dm

    if (stopifdecreasing & iter > 50) {
      if (E.log.lik < E.log.lik.old - 0.1) {
        #stop('Error: decreasing LL')
        continue  <- FALSE
        decreased <- TRUE
      }
    }
    if (iter/1000==round(iter/1000)) {
      cat('Iteration ',iter,' : ',Cm.diff,'  ',Dm.diff,'    ',E.log.lik,'    ',log.lik,'\n') # Cm.new.output$n.iter,'   '
    }
#    cat('\n\n')
#    cat(Cm.new.output$sigma2,'\t',Dm.new.output$sigma2)
#    cat('\n\n')
#    cat(Cm.new.output$W)
#    cat('\n\n')
#    cat(Dm.new.output$W)
#    cat('\n\n')

#    cat('\n\n')
#    cat('\n\n')
#    cat('\n\n')

    #############

    Cm <- Cm.new

    Dm <- Dm.new

    W.G <- W.G.new

    W.E <- W.E.new

    P.G <- P.G.new

    P.E <- P.E.new

    #############

    iter <- iter + 1

    #if (Cm.diff < tol & Dm.diff < tol) {continue <- F}
    if (abs(E.log.lik - E.log.lik.old) < tol.em) {continue <- F}

    # eigen(Cm)$values; eigen(Dm)$values
  }

  #################################

      if (compute.log.lik) {

      V.inv.array  <- make.V.inv.array(Vg=solve(Cm),Ve=solve(Dm),Dk=Dk)

      V.array      <- make.V.array(Vg=solve(Cm),Ve=solve(Dm),Dk=Dk)

      if (nc > 0) {X.transformed <- t(X) %*% Uk} else {X.transformed <- data.frame()}

      log.lik <- LL.diag(t(Y) %*% Uk,X=X.transformed,V.array=V.array,V.inv.array=V.inv.array)

      #H     <- kronecker(Dk,solve(Cm)) + kronecker(diag(n),solve(Dm))
      #H     <- kronecker(solve(Cm),K) + kronecker(diag(n),solve(Dm))

      # FULL:
      #H     <- kronecker(K,solve(Cm)) + kronecker(diag(n),solve(Dm))
      #Hinv  <- ginv(H)
      #bigX  <- kronecker(X,diag(p))
      #P     <- Hinv - Hinv %*% bigX %*% ginv(t(bigX) %*% Hinv %*% bigX) %*% t(bigX) %*% Hinv
      #log.lik <- LL(t(Y),H,P)

    } else {
      log.lik <- NA
    }

  #################################

  # prevent that Cm, Dm become asymmetric because of numerical inaccuracies
  Cm <- (Cm + t(Cm)) / 2
  Dm <- (Dm + t(Dm)) / 2

  if (is.null(rownames(Y))) {rownames(Y) <- paste0('genotype',1:n)}
  if (is.null(colnames(Y))) {colnames(Y) <- paste0('trait',1:p)}

  if (prediction & nc==0) {
    mu <- matrix(U %*% S.1 %*% t(Dm.sqrt.inv %*% Q.1),ncol=p)
  }

  pred.frame <- data.frame(trait=rep(colnames(Y),each=n),
                           genotype=rep(rownames(Y),p),
                           predicted=as.numeric(mu))

return(list(Cm=Cm,Dm=Dm,B=B,
            W.G=W.G,W.E=W.E,
            P.G=P.G,P.E=P.E,
            pred=pred.frame,log.lik=E.log.lik,log.lik2=log.lik,
            tol.em=tol.em,converged=(!continue),
            n=n,p=p,nc=nc,n.iter=iter,decreased=decreased))
}
#test <- EM_function(Y=Y,K=K,X=matrix(rep(1,nrow(K))),Dm.is.diagonal=F,pen.Cm=2,pen.Dm=2)

