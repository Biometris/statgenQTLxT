LL.grad.diag <- function(Y,X,Dk,V.inv.array) {

# X is the c x n covariate matrix, c being the number of covariates and
#      n being the number of genotypes
#      c has to be at least one (typically an intercept)
# Y is the p x n matrix of observed phenotypes, on p traits or environments
# x (optional) is an additional covariate.
#
# No missing values are allowed in any of X,Y and x
#
# It is assumed that X,Y and x have already been rotated by Uk, where Uk is such that
# the kinship matrix equals K = Uk %*% Dk %*% t(Uk)
# (R: w <- eigen(K); Dk <- diag(w$values); Uk <- w$vectors)
# Next, the original X, Y and x are post-multiplied by Uk, e.g. Y <- Y %*% Uk;
# see Zhou and Stephens 2014, supplement.
# It is the rotated versions that are the input of this function.


# Y=Yt;Dk=Dk;V.inv.array=V.inv.array

  stopifnot(ncol(Y)==ncol(X))

  nc      <- nrow(X)
  y       <- matrix(Y)
  p       <- nrow(Y)
  n       <- ncol(Y)

  q.scal <- sum(sapply(1:n,function(i){t(matrix(Y[,i])) %*% V.inv.array[i,,] %*% matrix(Y[,i])}))

  if (nc > 0) {
    q.vec  <- matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]),V.inv.array[i,,] %*% matrix(Y[,i]))})),1,sum))
    Q.matrix <-  matrix(apply((sapply(1:n,function(i){kronecker(matrix(X[,i]) %*% t(matrix(X[,i])),V.inv.array[i,,])})),1,sum),ncol=p)
  } else {
    q.vec  <- matrix(rep(0,p))
    Q.matrix <- matrix(rep(0,p*p),ncol=p)
  }

  #########################
  y <- matrix(Y)

  grad.G <- matrix(0,p,p)
  grad.E <- matrix(0,p,p)

  # Notation: see Zhou and Stephens 2014, supplement.
  #
  # We denote the quantities on the right hand side (RHS) of their equation (51) as follows:
  #
  # first term RHS (with italic q):  q.scal.G.i.j, q.scal.E.i.j
  #
  # second term RHS               : q.vec and Q.matrix, defined already above.
  #                                 the last factor we denote as q.vec.G.i.j, q.vec.E.i.j
  #
  # fourth term RHS               : we denote the matrix $Q_{ij}^g$ in the middle  as Q.G.i.j or Q.E.i.j
  #
  # GENERAL COMMENT: FIRST WE COMPUTE THE VARIANTS WITH 'E' (e.g. q.vec.E.i.j), and then for the corresponding
  #                  DERIVATIVE WITH RESPECT TO 'G' (e.g. q.vec.G.i.j) WE TYPICALLY ONLY HAVE TO MULTIPLY BY diag(Dk)

  Q.E.output <- array(dim=c(p,p,p,p))
  Q.G.output <- array(dim=c(p,p,p,p))

  Q.e.output <- array(dim=c(p,p,p,1))
  Q.g.output <- array(dim=c(p,p,p,1))

  for (j in 1:p) {
      for (i in 1:j) {
        ij.matrix2       <- make.ij.matrix(i=i,j=j,p=p,sym=F)
        Ind.i <- make.i.col.matrix(i=i,p=p)
        Ind.j <- make.i.col.matrix(i=j,p=p)

        #XX          <- sapply(1:n,function(m){matrix(X[,m]) %*% t(matrix(X[,m]))})# dim 1 is problematic !!!

        Q.E.i.j.matrix    <- sapply(1:n,function(m){kronecker(matrix(X[,m]) %*% t(matrix(X[,m])),matrix(V.inv.array[m,,i]) %*% t(matrix(V.inv.array[m,j,])))})
        Q.E.i.j           <- matrix(apply(Q.E.i.j.matrix,1,sum),ncol=p)
        Q.E.output[i,j,,] <-  Q.E.i.j

        q.vec.E.i.j.matrix <- sapply(1:n,function(m){kronecker(matrix(X[,m]),matrix(V.inv.array[m,,i]) %*% t(matrix(V.inv.array[m,j,])) %*% matrix(Y[,m]))})
        q.vec.E.i.j        <- matrix(apply(q.vec.E.i.j.matrix,1,sum))
        Q.e.output[i,j,,]  <-  q.vec.E.i.j

        q.scal.E.i.j.vec   <- sapply(1:n,function(m){t(matrix(Y[,m])) %*% matrix(V.inv.array[m,,i]) %*% t(matrix(V.inv.array[m,j,])) %*% matrix(Y[,m])})
        q.scal.E.i.j       <- sum(q.scal.E.i.j.vec)

        # Zhou and Stephens 2014, suppl, p.21, (42-44)
        Q.G.i.j     <- matrix(apply(t(t(Q.E.i.j.matrix) * diag(Dk)),1,sum),ncol=p)
        Q.G.output[i,j,,] <-  Q.G.i.j

        q.vec.G.i.j <- matrix(apply(t(t(q.vec.E.i.j.matrix) * diag(Dk)),1,sum))
        Q.g.output[i,j,,] <-  q.vec.G.i.j

        q.scal.G.i.j<- sum(q.scal.E.i.j.vec * diag(Dk))

        # Room for a speedup (?) : precompute x_{l} t(x_{l}), l=1...n and V_{l}^{-1} y_{l}
        # XX          <- sapply(1:n,function(m){matrix(X[,m]) %*% t(matrix(X[,m]))})# dim 1 is problematic !!!
        # Another speed-up could be achieved if the quantities Q.matrix, q.vec and q.scal
        # calculated at the beginning of this function are 'imported'; these are already calculated
        # in the functions LL.diag and LL.REML.diag

        # Zhou and Stephens 2014, suppl, p.21, (51)
        quad.form.G <- q.scal.G.i.j - t(q.vec) %*% solve(Q.matrix) %*% q.vec.G.i.j - t(q.vec.G.i.j) %*% solve(Q.matrix) %*% q.vec + t(q.vec) %*% solve(Q.matrix) %*% Q.G.i.j %*% solve(Q.matrix) %*% q.vec
        quad.form.E <- q.scal.E.i.j - t(q.vec) %*% solve(Q.matrix) %*% q.vec.E.i.j - t(q.vec.E.i.j) %*% solve(Q.matrix) %*% q.vec + t(q.vec) %*% solve(Q.matrix) %*% Q.E.i.j %*% solve(Q.matrix) %*% q.vec

        # Now compare with the 'inefficient' expression (32); there we have (I_{ij} + I_{ji}),
        # so the preceding terms need to be calculated again with i and j interchanged (if i!=j)
        if (i==j) {
          quad.form.G <- 2 * quad.form.G
          quad.form.E <- 2 * quad.form.E
        } else {
          Q.E.j.i.matrix <- sapply(1:n,function(m){kronecker(matrix(X[,m]) %*% t(matrix(X[,m])),matrix(V.inv.array[m,,j]) %*% t(matrix(V.inv.array[m,i,])))})
          Q.E.j.i     <- matrix(apply(Q.E.j.i.matrix,1,sum),ncol=p)
          Q.E.output[j,i,,] <-  Q.E.j.i

          #Q.E.j.i     <- t(Q.E.i.j)

          q.vec.E.j.i.matrix <- sapply(1:n,function(m){kronecker(matrix(X[,m]),matrix(V.inv.array[m,,j]) %*% t(matrix(V.inv.array[m,i,])) %*% matrix(Y[,m]))})
          q.vec.E.j.i <- matrix(apply(q.vec.E.j.i.matrix,1,sum))
          Q.e.output[j,i,,] <-  q.vec.E.j.i

          q.scal.E.j.i.vec <- sapply(1:n,function(m){t(matrix(Y[,m])) %*% matrix(V.inv.array[m,,j]) %*% t(matrix(V.inv.array[m,i,])) %*% matrix(Y[,m])})
          q.scal.E.j.i     <- sum(q.scal.E.j.i.vec) # q.scal.E.i.j

          # Zhou and Stephens 2014, suppl, p.21, (42-44)
          Q.G.j.i     <- matrix(apply(t(t(Q.E.j.i.matrix) * diag(Dk)),1,sum),ncol=p) #t(Q.G.i.j)
          Q.G.output[j,i,,] <-  Q.G.j.i

          q.vec.G.j.i <- matrix(apply(t(t(q.vec.E.j.i.matrix) * diag(Dk)),1,sum))
          Q.g.output[j,i,,] <-  q.vec.G.j.i

          q.scal.G.j.i<- sum(q.scal.E.j.i.vec * diag(Dk))

          quad.form.G <- quad.form.G + q.scal.G.j.i - t(q.vec) %*% solve(Q.matrix) %*% q.vec.G.j.i - t(q.vec.G.j.i) %*% solve(Q.matrix) %*% q.vec + t(q.vec) %*% solve(Q.matrix) %*% Q.G.j.i %*% solve(Q.matrix) %*% q.vec
          quad.form.E <- quad.form.E + q.scal.E.j.i - t(q.vec) %*% solve(Q.matrix) %*% q.vec.E.j.i - t(q.vec.E.j.i) %*% solve(Q.matrix) %*% q.vec + t(q.vec) %*% solve(Q.matrix) %*% Q.E.j.i %*% solve(Q.matrix) %*% q.vec

        }

        # Zhou and Stephens 2014, suppl, p.21, (54)
        trace.vec    <- sapply(1:n,function(m){(V.inv.array[m,i,j] + V.inv.array[m,j,i]) })
        trace.term.E <- sum(trace.vec)
        trace.term.G <- sum(trace.vec * diag(Dk))

        grad.G[i,j]      <- 0.5 * (quad.form.G - trace.term.G)
        grad.E[i,j]      <- 0.5 * (quad.form.E - trace.term.E)

      }
  }

  diag(grad.G) <- diag(grad.G) / 2
  diag(grad.E) <- diag(grad.E) / 2

return(list(grad=matrix(c(grad.G[upper.tri(grad.G,diag=T)],grad.E[upper.tri(grad.E,diag=T)])),
            Q.G.output=Q.G.output,Q.E.output=Q.E.output,
            Q.g.output=Q.g.output,Q.e.output=Q.e.output,
            q.vec=q.vec,Q.matrix=Q.matrix))
}
