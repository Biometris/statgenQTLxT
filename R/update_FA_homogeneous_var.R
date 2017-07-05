
update_FA_homogeneous_var <- function(Y=NULL,S=NULL,m,max.diag=1e04) {


  # Y : a matrix or data.frame dimensions n (observations) x p (variables)

  if (is.null(Y) & is.null(S)) {stop("Either the data(Y) or the sample covariance matrix(S) must be provided.")}

  if (!is.null(Y) & !is.null(S)) {stop("Either the data(Y) or the sample covariance matrix(S) must be provided.")}

  ################

  if (!is.null(Y)) {

    Y <- as.matrix(Y)

    if (sum(is.na(Y)) > 0) {stop('Y cannot contain missing values')}

    #col.means <- apply(Y,2,mean)

    Y <- scale(Y, center = TRUE, scale = FALSE)

    p <- ncol(Y)

    n <- nrow(Y)

    S <- t(Y) %*% Y / n

  } else {

    p <- ncol(S)

  }

  a <- eigen(S)

  ############

  if (m!=round(m)) {stop("m needs to be a positive integer")}
  if (m < 1) {stop("m needs to be a positive integer")}
  if (m >= p) {stop("m needs to be smaller than the number of variables")}

  #########

  sigma2 <- mean(a$values[-(1:m)])

  if (sigma2 < 1/max.diag) {sigma2 <- 1/max.diag}

  if (m == 1) {
    W <- matrix(a$vectors[,1:m]) %*% matrix(sqrt(a$values[1:m] - sigma2))
  } else {
    W <- a$vectors[,1:m] %*% diag(sqrt(a$values[1:m] - sigma2))
  }
return(list(W=W,sigma2=sigma2))
}

#

