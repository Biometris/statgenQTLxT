update_LS <- function(Y,rho,
                      W.start=NULL,
                      m=ifelse(is.null(W.start), yes=2, no=ncol(W.start)),
                      P.start=NULL,
                      tol.em=0.0001,
                      max.iter=100,
                      print.progress=FALSE,
                      penalize.diagonal=TRUE) {

  # to do : extension to Y- mu; now mu is assumed to be zero

  # Y=t(cObs);rho=10;W.start=NULL;P.start=NULL;m=2;tol.em=0.001;max.iter=1000;print.progress=T;penalize.diagonal=F
  #  Y=A; rho=pen.Cm; W.start=W.G; P.start=P.G; penalize.diagonal=penalize.diagonal.G ; m=ncol(W.start); print.progress=T

  # Y : a matrix or data.frame dimensions n (observations) x p (variables)

  ################

  Y <- as.matrix(Y)

  if (sum(is.na(Y)) > 0) {stop('Y cannot contain missing values')}

  col.means <- apply(Y,2,mean)

  Y <- scale(Y, center = TRUE, scale = FALSE)

  p <- ncol(Y)

  n <- nrow(Y)

  ############

  if (!is.null(P.start)) {
    stopifnot(class(P.start) %in% c('matrix','data.frame'))
    stopifnot(ncol(P.start)==p & nrow(P.start)==p)
  }

  if (!is.null(W.start)) {
    stopifnot(class(W.start) %in% c('matrix','data.frame'))
    if (ncol(W.start)!=m) {stop('m needs to be equal to the number of columns of W.start')}
    if (is.null(P.start)) {stop('W.start and P.start should be either both NULL (default), or both have a sensible value')}
    stopifnot(nrow(W.start)==p)
  } else {
    if (!is.null(P.start)) {stop('W.start and P.start should be either both NULL (default), or both have a sensible value')}
  }

  if (m!=round(m)) {stop("m needs to be integer")}
  if (m < 0) {stop("m cannot be negative")}
  if (m >= p) {stop("m needs to be smaller than the number of variables")}


  ############

  if (is.null(W.start)) {

    S <- t(Y) %*% Y / n

    a <- eigen(S)

    sigma2 <- mean(a$values[-(1:m)])

    P.start <- diag(p) / sigma2

    W.start <- a$vectors[,1:m] %*% diag(sqrt(a$values[1:m] - sigma2))

  }

  ############

  W   <- W.start

  P   <- P.start

  ############

  continue <- T

  iter <- 1

  total.diff <- 10000000000

  #############################################
  # EM

  while (continue & iter < max.iter) {

    # prevent that P become asymmetric because of numerical inaccuracies
    P <- (P + t(P)) / 2        # p x p

    B <- t(W) %*% P %*% W      # m x m

    Sigma <- ginv(diag(m) + B) # m x m

    M1 <-  Sigma %*% t(W) %*% P %*% t(Y) # m x n

    A <- ginv(n * Sigma + M1 %*% t(M1))  # m x m

    W.new <- t(Y) %*% t(M1) %*% A  # p x m

    new.data <- t(Y) - W.new %*% M1 # p x n

    new.S <- new.data %*% t(new.data) / n

    P.new <- glasso(s=new.S, rho=rho,penalize.diagonal=penalize.diagonal)$wi

    P.diff <- sum(abs(P.new - P))

    W.diff <- sum(abs(W.new - W))

    total.diff <- P.diff + W.diff

    if (print.progress) {
      cat('Iteration ',iter,' : ',P.diff,'  ',W.diff,'\n')
    }

    P <- P.new

    W <- W.new

    iter <- iter + 1

    if (total.diff < tol.em) {continue <- F}

  }
#image(P)

return(list(W=W,P=P,n.iter=iter))
}

