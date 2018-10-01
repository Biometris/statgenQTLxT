#include <RcppArmadillo.h>
#include "getThr.h"

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
List fastGLSIBDCPP(const arma::cube &mp,
                   const arma::vec &y,
                   const arma::mat &sigma,
                   unsigned int ref,
                   Rcpp::Nullable<Rcpp::NumericVector> size_param = R_NilValue,
                   Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  // Check that reference allele contains not only zeros.
  if (all(vectorise(mp.slice(ref - 1)) == 0)) {
    stop("Invalid reference allele.\n"
         "Assure the reference allele contains non-zero values.\n");
  }
  // Get number of genotypes, alleles and markers.
  unsigned int m = mp.n_slices - 1;
  unsigned int p = mp.n_cols;
  unsigned int n = mp.n_rows;
  // Define covs as intercept.
  arma::mat covs = ones<mat>(n, 1);
  if (size_param.isNotNull()) {
    // Add other covariates.
    covs.insert_cols(1, Rcpp::as<arma::mat>(size_param));
  }
  // Get number of covariates (including intercept).
  unsigned int nCov = covs.n_cols;
  // Compute inverse of t(M).
  arma::mat mt = chol(sigma, "lower").i();
  // pre-multiply the IBD probablities with t(M).
  // Skip reference allele and move alleles past reference allele to left.
  arma::cube mp2 = cube(n, p, m);
  int nThr = getThr(nCores);
#pragma omp parallel for num_threads(nThr)
  for (unsigned int i = 0; i < m + 1; i++) {
    if (i < ref - 1) {
      mp2.slice(i) = mt * mp.slice(i);
    } else if (i > ref - 1) {
      mp2.slice(i - 1) = mt * mp.slice(i);
    }
  }
  // Pre-multiply the phenotype (y) with t(M).
  arma::mat tMy = mt * y;
  // pre-multiply the intercept and covariates with t(M).
  arma::mat tMfixCovs = mt * covs;
  // Compute residuals and RSS over all markers.
  arma::mat Q, R;
  arma::qr_econ(Q, R, tMfixCovs);
  double RSSEnv = accu(square(tMy - Q * Q.t() * tMy));
  // Define matrices for storing output.
  arma::mat beta1 = mat(tMfixCovs.n_cols, p);
  arma::mat beta2 = zeros<mat>(m, p);
  NumericVector pVal = NumericVector(p);
  NumericVector RLR2 = NumericVector(p);
  // Loop over markers and compute betas and RSS.
#pragma omp parallel for num_threads(nThr)
  for (unsigned int i = 0; i < p; i ++) {
    arma::mat X = mp2(span(), span(i), span());
    // Get indices of alleles in X that are not entirely 0.
    arma::uvec posInd = find( sum( abs(X), 0 ) > 0 );
    // Subset X on those columns
    X = X.cols(posInd);
    // Algorithm for computing beta1, beta2 and RSSFull for current marker.
    arma::mat tX2VinvX2Inv = (X.t() * X).i();
    arma::mat tX1VinvX2 = tMfixCovs.t() * X;
    arma::mat tX2VinvY = X.t() * tMy;
    arma::mat XS = tMfixCovs.t() * tMfixCovs -
      tX1VinvX2 * tX2VinvX2Inv * tX1VinvX2.t();
    arma::vec beta1j = solve(XS, tMfixCovs.t() * tMy -
      tX1VinvX2 * tX2VinvX2Inv * tX2VinvY);
    arma::vec beta2j = tX2VinvX2Inv * (tX2VinvY - tX1VinvX2.t() * beta1j);
    double RSSFullj = accu(square(tMy - tMfixCovs * beta1j - X * beta2j));
    // Put results for beta1 for current marker in output matrix.
    beta1.col(i) = beta1j;
    // Create vector of zeros for beta2. Not possible before computing
    // RSSFullj since dimensions won't match for X and beta2j then.
    arma::vec beta2jZeros = zeros<vec>(m);
    // Fill non-zero position with computed values and add to output matrix.
    beta2jZeros(posInd) = beta2j;
    beta2.col(i) = beta2jZeros;
    // Compute further output elements.
    double df1x = beta2j.n_elem;
    double df2x = n - df1x - nCov;
    double fValx = (RSSEnv - RSSFullj) / RSSFullj * df2x / df1x;
    pVal(i) = R::pf(fValx, df1x, df2x, false, false);
    // Compute R_LR^2 statistic from Sun et al 2010, heredity.
    RLR2[i] = 1 - exp((RSSFullj - RSSEnv) / n);
  }
  // Create output.
  return List::create(_["beta1"] = beta1,
                      _["beta2"] = beta2,
                      _["pVal"] = pVal,
                      _["RLR2"] = RLR2);
}
