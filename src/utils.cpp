#include <RcppArmadillo.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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
  eigVecs = fliplr( eigVecs.tail_cols(nPca) );
  arma::mat S = arma::diagmat( reverse(eigVals.tail(nPca)) );
  return eigVecs * S * eigVecs.t();
}
