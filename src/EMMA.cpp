#include <RcppArmadillo.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' EMMA helper functions
//'
//' Helper functions for computing REML estimates of genetic and residual
//' variance components using the EMMA algorithm.
//'
//' @inheritParams EMMA
//' @param X a q x n covariate matrix, q being the number of covariates and n
//' being the number of genotypes. q has to be at least one (typically an
//' intercept).
//'
//' @keywords internal
// [[Rcpp::export]]
List emmaEigenR(const arma::mat k,
                const arma::mat x) {
  int n = x.n_rows;
  int q = x.n_cols;
  // Compute n-q non-zero eigenvalues of SHS as defined in eqn. 5 of Kang.
  arma::mat s = arma::eye(n, n);
  if (q == 1) {
    s.for_each( [ n ](mat::elem_type& val) { val -= 1 / double (n); } );
  } else {
    s -= x * solve(x.t() * x, x.t());
  }
  arma::vec eigVals(n);
  arma::mat eigVecs(n, n);
  arma::eig_sym(eigVals, eigVecs, s * (k + eye(n, n)) * s);
  return List::create(_["values"] = reverse(eigVals.tail(n - q)) - 1,
                      _["vectors"] = fliplr( eigVecs.tail_cols(n - q) ));
}
