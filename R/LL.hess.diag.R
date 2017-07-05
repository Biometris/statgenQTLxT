LL.hess.diag <- function(Y,X,Dk,V.inv.array,Q.G.output,Q.E.output,
                                            Q.g.output,Q.e.output,
                                            q.vec,Q.matrix) {

  #Y=Yt

  dk <- diag(Dk)

  dk2 <- dk^2

  p <- nrow(Y)

  n       <- ncol(Y)

  Q <- p * (p + 1) / 2

  y <- matrix(Y)

  hess.G <- matrix(0,Q,Q)

  hess.E <- matrix(0,Q,Q)

  hess.GE <- matrix(0,Q,Q)

  for (g2 in 1:Q) {
    for (g1 in 1:g2) {

      # g1 <- 1; g2 <- 2

      g1.indices <- which.i(i=g1,p=Q)
      g2.indices <- which.i(i=g2,p=Q)

      i          <- g1.indices[1]
      j          <- g1.indices[2]

      ip         <- g2.indices[1]
      jp         <- g2.indices[2]

      prod.mat.T <- sapply(1:n,function(m){(as.numeric(V.inv.array[m,j,ip] * matrix(V.inv.array[m,,i]) %*% t(matrix(V.inv.array[m,jp,])))+as.numeric(V.inv.array[m,j,jp] * matrix(V.inv.array[m,,i]) %*% t(matrix(V.inv.array[m,ip,])))+as.numeric(V.inv.array[m,i,ip] * matrix(V.inv.array[m,,j]) %*% t(matrix(V.inv.array[m,jp,])))+as.numeric(V.inv.array[m,i,jp] * matrix(V.inv.array[m,,j]) %*% t(matrix(V.inv.array[m,ip,])))) })
      prod.mat.T <- array(t(prod.mat.T),dim=c(n,p,p))

      Qinv       <- solve(Q.matrix)

      ##############################################

      q.scal.EE.i.j.ip.jp.vec   <- sapply(1:n,function(m){t(matrix(Y[,m])) %*% prod.mat.T[m,,] %*% matrix(Y[,m])})
      q.scal.EE.i.j.ip.jp       <- sum(q.scal.EE.i.j.ip.jp.vec)

      Q.EE.i.j.ip.jp.matrix     <- sapply(1:n,function(m){kronecker(matrix(X[,m]) %*% t(matrix(X[,m])),prod.mat.T[m,,])})
      Q.EE.i.j.ip.jp            <- matrix(apply(Q.EE.i.j.ip.jp.matrix,1,sum),ncol=p)

      q.vec.EE.i.j.ip.jp.matrix <- sapply(1:n,function(m){kronecker(matrix(X[,m]),prod.mat.T[m,,] %*% matrix(Y[,m]))})
      q.vec.EE.i.j.ip.jp        <- matrix(apply(q.vec.EE.i.j.ip.jp.matrix,1,sum))

      q.vec.EE.jp.ip.j.i.matrix <- sapply(1:n,function(m){V.inv.array[m,ip,j] * kronecker(matrix(X[,m]),matrix(V.inv.array[m,,jp]) %*% t(matrix(V.inv.array[m,i,])) %*% matrix(Y[,m])) + V.inv.array[m,jp,j] * kronecker(matrix(X[,m]),matrix(V.inv.array[m,,ip]) %*% t(matrix(V.inv.array[m,i,])) %*% matrix(Y[,m])) + V.inv.array[m,ip,i] * kronecker(matrix(X[,m]),matrix(V.inv.array[m,,jp]) %*% t(matrix(V.inv.array[m,j,])) %*% matrix(Y[,m])) + V.inv.array[m,jp,i] * kronecker(matrix(X[,m]),matrix(V.inv.array[m,,ip]) %*% t(matrix(V.inv.array[m,j,])) %*% matrix(Y[,m]))})
      q.vec.EE.jp.ip.j.i        <- matrix(apply(q.vec.EE.jp.ip.j.i.matrix,1,sum))

      ####################################

      q.scal.GG.i.j.ip.jp       <- sum(dk2 * q.scal.EE.i.j.ip.jp.vec)

      Q.GG.i.j.ip.jp            <- matrix(apply(t(dk2 * t(Q.EE.i.j.ip.jp.matrix)),1,sum),ncol=p)

      q.vec.GG.i.j.ip.jp        <- matrix(apply(t(dk2 * t(q.vec.EE.i.j.ip.jp.matrix)),1,sum))

      q.vec.GG.jp.ip.j.i        <- matrix(apply(t(dk2 * t(q.vec.EE.jp.ip.j.i.matrix)),1,sum))

      ####################################

      q.scal.GE.i.j.ip.jp       <- sum(dk * q.scal.EE.i.j.ip.jp.vec)

      Q.GE.i.j.ip.jp            <- matrix(apply(t(dk * t(Q.EE.i.j.ip.jp.matrix)),1,sum),ncol=p)

      q.vec.GE.i.j.ip.jp        <- matrix(apply(t(dk * t(q.vec.EE.i.j.ip.jp.matrix)),1,sum))

      q.vec.GE.jp.ip.j.i        <- matrix(apply(t(dk * t(q.vec.EE.jp.ip.j.i.matrix)),1,sum))

      ####################################

      t1 <- q.scal.EE.i.j.ip.jp

      t2 <-  -1 * as.numeric(t(q.vec) %*% Qinv %*% q.vec.EE.i.j.ip.jp)

      t3 <-  -1 * as.numeric(t(q.vec.EE.jp.ip.j.i) %*% Qinv %*% q.vec)

      t4 <-  -1 * as.numeric(t(matrix(Q.e.output[j,i,,]) + matrix(Q.e.output[i,j,,])) %*% Qinv %*% (matrix(Q.e.output[ip,jp,,]) + matrix(Q.e.output[jp,ip,,])))

      t5 <- as.numeric(t(q.vec) %*% Qinv %*% (Q.E.output[i,j,,] + Q.E.output[j,i,,]) %*% Qinv %*% (matrix(Q.e.output[ip,jp,,]) + matrix(Q.e.output[jp,ip,,])))

      t6 <- as.numeric(t(matrix(Q.e.output[j,i,,]) + matrix(Q.e.output[i,j,,]))) %*% Qinv %*% (Q.E.output[ip,jp,,] + Q.E.output[jp,ip,,]) %*% Qinv %*% q.vec

      t7 <- as.numeric(t(q.vec) %*% Qinv %*% Q.EE.i.j.ip.jp  %*% Qinv %*% q.vec)

      t8 <- -1 * as.numeric(t(q.vec) %*% Qinv %*% (Q.E.output[i,j,,] + Q.E.output[j,i,,]) %*% Qinv %*% (Q.E.output[ip,jp,,] + Q.E.output[jp,ip,,]) %*% Qinv %*% q.vec)

      ############################################

      t1.GG <- q.scal.GG.i.j.ip.jp

      t2.GG <-  -1 * as.numeric(t(q.vec) %*% Qinv %*% q.vec.GG.i.j.ip.jp)

      t3.GG <-  -1 * as.numeric(t(q.vec.GG.jp.ip.j.i) %*% Qinv %*% q.vec)

      t4.GG <-  -1 * as.numeric(t(matrix(Q.g.output[j,i,,]) + matrix(Q.g.output[i,j,,])) %*% Qinv %*% (matrix(Q.g.output[ip,jp,,]) + matrix(Q.g.output[jp,ip,,])))

      t5.GG <- as.numeric(t(q.vec) %*% Qinv %*% (Q.G.output[i,j,,] + Q.G.output[j,i,,]) %*% Qinv %*% (matrix(Q.g.output[ip,jp,,]) + matrix(Q.g.output[jp,ip,,])))

      t6.GG <- as.numeric(t(matrix(Q.g.output[j,i,,]) + matrix(Q.g.output[i,j,,]))) %*% Qinv %*% (Q.G.output[ip,jp,,] + Q.G.output[jp,ip,,]) %*% Qinv %*% q.vec

      t7.GG <- as.numeric(t(q.vec) %*% Qinv %*% Q.GG.i.j.ip.jp  %*% Qinv %*% q.vec)

      t8.GG <- -1 * as.numeric(t(q.vec) %*% Qinv %*% (Q.G.output[i,j,,] + Q.G.output[j,i,,]) %*% Qinv %*% (Q.G.output[ip,jp,,] + Q.G.output[jp,ip,,]) %*% Qinv %*% q.vec)

      ############################################

      t1.GE <- q.scal.GE.i.j.ip.jp

      t2.GE <-  -1 * as.numeric(t(q.vec) %*% Qinv %*% q.vec.GE.i.j.ip.jp)

      t3.GE <-  -1 * as.numeric(t(q.vec.GE.jp.ip.j.i) %*% Qinv %*% q.vec)

      t7.GE <- as.numeric(t(q.vec) %*% Qinv %*% Q.GE.i.j.ip.jp  %*% Qinv %*% q.vec)

      t4.GE <-  -1 * as.numeric(t(matrix(Q.g.output[j,i,,]) + matrix(Q.g.output[i,j,,])) %*% Qinv %*% (matrix(Q.e.output[ip,jp,,]) + matrix(Q.e.output[jp,ip,,])))

      t5.GE <- as.numeric(t(q.vec) %*% Qinv %*% (Q.G.output[i,j,,] + Q.G.output[j,i,,]) %*% Qinv %*% (matrix(Q.e.output[ip,jp,,]) + matrix(Q.e.output[jp,ip,,])))

      t6.GE <- as.numeric(t(matrix(Q.g.output[j,i,,]) + matrix(Q.g.output[i,j,,]))) %*% Qinv %*% (Q.E.output[ip,jp,,] + Q.E.output[jp,ip,,]) %*% Qinv %*% q.vec

      t8.GE <- -1 * as.numeric(t(q.vec) %*% Qinv %*% (Q.G.output[i,j,,] + Q.G.output[j,i,,]) %*% Qinv %*% (Q.E.output[ip,jp,,] + Q.E.output[jp,ip,,]) %*% Qinv %*% q.vec)

      ############################################

      E.quad.form        <- t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8

      E.trace            <- sum(V.inv.array[,jp,i] * V.inv.array[,j,ip]) + sum(V.inv.array[,ip,i] * V.inv.array[,j,jp]) + sum(V.inv.array[,jp,j] * V.inv.array[,i,ip]) + sum(V.inv.array[,ip,j] * V.inv.array[,i,jp])

      G.trace            <- sum(dk2 * V.inv.array[,jp,i] * V.inv.array[,j,ip]) + sum(dk2 * V.inv.array[,ip,i] * V.inv.array[,j,jp]) + sum(dk2 * V.inv.array[,jp,j] * V.inv.array[,i,ip]) + sum(dk2 * V.inv.array[,ip,j] * V.inv.array[,i,jp])

      G.quad.form        <- t1.GG + t2.GG + t3.GG + t4.GG + t5.GG + t6.GG + t7.GG + t8.GG

      GE.trace            <- sum(dk * V.inv.array[,jp,i] * V.inv.array[,j,ip]) + sum(dk * V.inv.array[,ip,i] * V.inv.array[,j,jp]) + sum(dk * V.inv.array[,jp,j] * V.inv.array[,i,ip]) + sum(dk * V.inv.array[,ip,j] * V.inv.array[,i,jp])

      GE.quad.form        <- t1.GE + t2.GE + t3.GE + t4.GE + t5.GE + t6.GE + t7.GE + t8.GE

      hess.E[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * E.trace - E.quad.form)

      hess.G[g1,g2]      <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * G.trace - G.quad.form)

      hess.GE[g1,g2]     <- (1 / ((1 + as.numeric(i==j))*(1 + as.numeric(ip==jp)))) * (0.5 * GE.trace - GE.quad.form)

    }
  }

  hess.G <- hess.G + t(hess.G); diag(hess.G) <- diag(hess.G) / 2
  hess.E <- hess.E + t(hess.E); diag(hess.E) <- diag(hess.E) / 2
  hess.GE<- hess.GE+ t(hess.GE);diag(hess.GE)<- diag(hess.GE)/ 2

  hessian <- rbind(cbind(hess.G,hess.GE),cbind(t(hess.GE),hess.E))

return(hessian)
}
