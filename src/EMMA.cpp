#include <RcppArmadillo.h>

#define _USE_MATH_DEFINES // for C++
#include <cmath>

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
  eigVals = reverse(eigVals.tail(n - q)) - 1;
  return List::create(_["values"] = arma::conv_to< std::vector<double>>::from(eigVals),
                      _["vectors"] = fliplr( eigVecs.tail_cols(n - q) ));
}

// [[Rcpp::export]]
std::vector<double> emmaREMLLL(double logDelta,
                               arma::vec lambda,
                               arma::vec etas1,
                               double n,
                               double t,
                               arma::vec etas2) {
  // Compute the REML LL as in eqn. 7 of Kang.
  double nq = etas1.n_elem + n - t;
  double delta = exp(logDelta);
  lambda += delta;
  arma::vec ll = 0.5 * (nq * (log(nq / (2 * M_PI)) - 1 - log(sum(square(etas1) /
    lambda) + etas2 / delta)) - sum(log(lambda)) + (t - n) * logDelta);
  return arma::conv_to< std::vector<double>>::from(ll);
}
