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
arma::mat multiAllKin(arma::cube x,
                      arma::vec posCor,
                      Rcpp::Nullable<Rcpp::NumericVector> denom = R_NilValue) {
  // Construct K by summing tcrossprod over alleles per marker.
  // Multiply by correction for position.
  arma::mat k = zeros<mat>(x.n_rows, x.n_rows);
  for (unsigned int m = 0; m < x.n_cols; m++) {
    arma::mat xm = x(span(), span(m), span());
    k += xm * xm.t() * posCor(m);
  }
  // To get ones on diagonal of final matrix put sum of position corrections.
  k.diag() = sum(posCor) * ones<vec>(x.n_rows);
  // Compute denominator.
  double denominator;
  if (denom.isNull()) {
    denominator = sum(posCor);
  } else {
    denominator = Rcpp::as<double>(denom);
  }
  // Divide by sum of position correction - almost the chr length exept for the
  // extra bits for first and last marker.
  return k / denominator;
}
