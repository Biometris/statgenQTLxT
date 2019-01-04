#include <RcppArmadillo.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Compute square root of a symmetric, positive definite matrix
//'
//' Given a symmetric, positive definite matrix X a matrix Y is computed such
//' that \eqn{Y^2 = X}. Computation is done using eigendecomposition of X.
//'
//' @param X A symmetric, positive definite matrix.
//'
//' @return A matrix Y such that \eqn{Y^2 = X}.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat matrixRoot(const arma::mat x) {
  try {
    return sqrtmat_sympd(x);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error("x should be a symmetric positive definite matrix.\n");
  }
}

//' Reduce the kinship matrix
//'
//' The kinship matrix is reduced using nPca eigenvectors of K.
//'
//' @inheritParams runMultiTraitGwas
//'
//' @param nPca An integer, the number of eigenvectors used for reducing the
//' kinship matrix.
//'
//' @return The reduced kinship matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat reduceKinship(const arma::mat K,
                        const int nPca) {
  arma::vec eigVals( K.n_cols );
  arma::mat eigVecs( size(K) );
  arma::eig_sym(eigVals, eigVecs, K);
  // fliplr and reverse are needed because eigenVals and eigenVecs are returned
  // sorted in ascending order.
  eigVecs = fliplr( eigVecs.tail_cols(nPca) );
  arma::mat S = arma::diagmat( reverse(eigVals.tail(nPca)) );
  return eigVecs * S * eigVecs.t();
}
