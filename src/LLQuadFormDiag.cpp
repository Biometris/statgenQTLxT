#include <RcppArmadillo.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double LLQuadFormDiagCPP(const arma::mat &Y,
                         const arma::cube &vInv,
                         Rcpp::Nullable<Rcpp::NumericVector> size_param = R_NilValue) {
  unsigned int nc = 0;
  arma::mat X;
  // Define covs as intercept.
  if (size_param.isNotNull()) {
    // Add other covariates.
    X = Rcpp::as<arma::mat>(size_param);
    nc = X.n_rows;
  }
  // Get number of genotypes and markers.
  unsigned int n = Y.n_cols;
  unsigned int p = Y.n_rows;
  // Compute scalair part.
  double qScal = 0;
  for (unsigned int i = 0; i < n; i++) {
    qScal += as_scalar(Y.col(i).t() * vInv.slice(i) * Y.col(i));
  }
  if (nc > 0) {
    arma::vec qVec = zeros<vec>(p);
    arma::mat qMat = zeros<mat>(p, p * nc);
    // compute q and Q.
    for (unsigned int i = 0; i < n; i++) {
      qVec += kron(X.col(i), vInv.slice(i) * Y.col(i));
      qMat += kron(X.col(i) * X.col(i).t(), vInv.slice(i));
    }
    // Compute quadratic part.
    arma::mat qQuad = qVec.t() * solve(qMat, qVec);
    qScal -= as_scalar(qQuad * qQuad.t());
  }
  return qScal;
}


// LLQuadFormDiag <- function(Y,
//                            X = data.frame(),
//                            vInvLst) {
//   nc <- nrow(X)
//   n <- ncol(Y)
//   p <- nrow(Y)
// ## Define function for faster computation of scalar part.
//   scalFunc <- function(i) {
//     as.numeric(Matrix::crossprod(Y[, i, drop = FALSE],
//                                  vInvLst[[i]] %*% Y[, i, drop = FALSE]))
//   }
// ## Compute scalair part.
//   qScal <- sum(sapply(X = 1:n, FUN = scalFunc))
//   if (nc > 0) {
// ## Define functions for faster computation of q and Q.
//     qVecFunc <- function(i) {
//       as.numeric(Matrix::kronecker(X[, i, drop = FALSE],
//                                    vInvLst[[i]] %*% Y[, i, drop = FALSE]))
//     }
//     qMatFunc <- function(i) {
//       as.numeric(Matrix::kronecker(Matrix::tcrossprod(X[, i, drop = FALSE]),
//                                    vInvLst[[i]]))
//     }
// ## Compute q and Q.
//     if (p == 1 && nc == 1) {
//       qVec <- sum(sapply(X = 1:n, FUN = qVecFunc))
//       QMatrix <- sum(sapply(X = 1:n, FUN = qMatFunc))
//     } else {
//       qVec  <- rowSums(sapply(X = 1:n, FUN = qVecFunc))
//       QMatrix <- matrix(rowSums(sapply(X = 1:n, FUN = qMatFunc)), ncol = p * nc)
//     }
// ## Compute quadratic part.
//     quadFormPart <- qScal - as.numeric(crossprod(qVec %*% solve(QMatrix, qVec)))
//   } else {
//     quadFormPart <- qScal
//   }
//   return(quadFormPart)
// }

