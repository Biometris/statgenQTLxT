#' Fast Generalized Least Squares algoritm for IBD based QTL mapping
#'
#' Compute statistics for the Generalized Least Squares (GLS) F-test for
#' IBD based QTL Mapping.
#'
#' @inheritParams fastGLS
#'
#' @param MP An array of marker probabilities.
#' @param ref An integer indicating the allele to use as reference allele in
#' the computations.
#'
#' @keywords internal
fastGLSIBD <- function(y,
                       MP,
                       Sigma,
                       covs,
                       ref) {
  ## Check class and missing values.
  if (missing(y) || !(inherits(y, "Matrix") || is.numeric(y)) || anyNA(y)) {
    stop("y should be a numeric vector without missing values.\n")
  }
  if (missing(Sigma) || !(inherits(Sigma, "Matrix") || is.matrix(Sigma)) ||
      anyNA(Sigma)) {
    stop("Sigma should be a matrix without missing values.\n")
  }
  if (!is.null(covs) && (!(inherits(covs, "Matrix") || is.matrix(covs)) ||
                         anyNA(covs))) {
    stop("covs should be a numeric vector without missing values.\n")
  }
  n <- length(y)
  ## Check dimensions.
  if (nrow(Sigma) != n || ncol(Sigma) != n) {
    stop(paste("The number of elements in y should be identical to the",
               "number of rows and columns in Sigma.\n"))
  }
  if (!is.null(covs) && nrow(covs) != n) {
    stop(paste("The number of elements in y should be identical to the",
               "number of rows in covs.\n"))
  }
  resCpp <- fastGLSIBDCPP(MP, y, Sigma, ref, covs, ncores = 4)
  beta1 <- resCpp$beta1
  beta2 <- resCpp$beta2
  FVal <- resCpp$FVal
  df1 <- resCpp$df1
  df2 <- resCpp$df2
  pVal <- pf(q = FVal, df1 = df1, df2 = df2, lower.tail = FALSE)
  return(list(beta1 = beta1, beta2 = beta2, pvalues = pVal))
}

#' fastGLS IBD Original function
#'
#' @keywords internal
fastGLSIBDOrig <- function(V,
                           Xcov,
                           y,
                           MP,
                           ref,
                           ind.matrix) {
  m <- dim(MP)[3]
  p <- dim(MP)[2]
  n <- length(y)
  nc <- ncol(Xcov)
  Vinv <- solve(V)

  #############################################################
  # just to test the 1st marker; only SSred is still used later on
  #X <- cbind(Xcov, MP[,1, -(ref)])

  # vector of GLS estimates, for the intercept and the non-reference alleles
  #b <- solve(t(X)%*% Vinv %*% X) %*% t(X)%*% Vinv %*% y

  # compute sums of squares, for full and reduced model
  #P <- Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv
  Pr <- Vinv - Vinv %*% Xcov %*% solve(t(Xcov) %*% Vinv %*% Xcov) %*% t(Xcov) %*% Vinv
  #SSfull <- t(y) %*% P %*% y

  # Also Used below:
  SSred <- t(y) %*% Pr %*% y

  # up to a constant ...
  #Fstat <- (SSred - SSfull) / SSfull

  ############################################################
  # fast version ...

  date()

  Fstat_vector <- rep(NA, p)

  pvalues      <- rep(NA, p)

  M <- Matrix::solve(Matrix::chol(V))
  #max(abs(Vinv - M%*%t(M)))

  ## Pre-multiply the phenotype (y) with t(M).
  tMy <- Matrix::crossprod(M, y)
  ## pre-multiply the intercept and covariates with t(M).
  tMfixCovs <- Matrix::crossprod(M, Xcov)

  # reduced form, without the ref allele
  tMPr  <- MP[,,-(ref)]

  ## pre-multiply the IBD probablities with t(M)
  for (k in 1:(m - 1)) {
    tMPr[,,k] <- Matrix::crossprod(M, tMPr[,,k])
  }

  # The idea is now to apply the standard ordinary least squares
  # formulas for a partitioned design matrix X = [X1 X2], where X1 = tMfixCovs
  # are the transformed fixed covariates and X2 = one of the loci in tMPr
  # (the transformed probabilities)

  #max(abs(tMPr - MP[,,-(ref)]))

  # Denoting the design matrix for a single locus as Mj,
  # compute Mj^t Mj, for all loci
  tMPr2     <- lapply(X = 1:p, FUN = function(j) tMPr[,j , ])
  tX2VinvX2 <- lapply(tMPr2, crossprod)

  tX1VinvX2 <- lapply(tMPr2, crossprod, x = tMfixCovs)
  tX2VinvY  <- lapply(tMPr2, crossprod, y = tMy)

  tX1VinvX1 <- t(tMfixCovs) %*% tMfixCovs
  tX1VinvY  <- t(tMfixCovs) %*% tMy

  # store the results in these matrices, corresponding to the
  # fixed covariates,  and the m-1 alleles
  beta1 <- matrix(0,  nc, p)
  beta2 <- matrix(0, m - 1, p)

  # perform QTL-mapping

  for (j in 1:p) {
    #j=1
    nm <- which(ind.matrix[,j] == 1)

    tX2VinvX2_inv_j <- solve(tX2VinvX2[[j]][nm, nm, drop = FALSE])

    tX1VinvX2_j <- tX1VinvX2[[j]][ , nm, drop = FALSE]

    tX2VinvY_j <- tX2VinvY[[j]][nm, , drop = FALSE]


    X_S_inv_j <- solve(tX1VinvX1 - (tX1VinvX2_j %*% tX2VinvX2_inv_j %*% t(tX1VinvX2_j)))

    beta1[,j] <- as.numeric(X_S_inv_j %*% tX1VinvY - X_S_inv_j %*% tX1VinvX2_j %*% tX2VinvX2_inv_j %*% tX2VinvY_j)


    beta2[nm, j] <- as.numeric(tX2VinvX2_inv_j %*% tX2VinvY_j + tX2VinvX2_inv_j %*% t(tX1VinvX2_j) %*% X_S_inv_j %*% tX1VinvX2_j %*% tX2VinvX2_inv_j %*% tX2VinvY_j - tX2VinvX2_inv_j %*% t(tX1VinvX2_j) %*% X_S_inv_j %*% tX1VinvY)

    residual_j <- as.numeric(tMy - tMfixCovs %*% as.matrix(beta1[,j]) - tMPr2[[j]] %*% as.matrix(beta2[,j]))

    SSfull_j <- sum(residual_j * residual_j)

    df1 <- length(nm)
    df2 <- n - length(nm) - nc

    # SSred is obtained from the 'undiagonalized' code at the beginning
    Fstat_vector[j] <- (df2 / df1) * (SSred - SSfull_j) / SSfull_j

    #FVal <- (RSSEnv - RSSFull) / RSSFull * df2
    pvalues[j] <- pf(q = Fstat_vector[j], df1 = df1, df2 = df2, lower.tail = FALSE)


  }
  return(list(beta1 = beta1, beta2 = beta2, pvalues = pvalues))
}




