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
List fastGLSIBDCPP(arma::cube MP,
                   arma::vec y,
                   arma::mat sigma,
                   int ref,
                   Rcpp::Nullable<Rcpp::NumericVector> size_param = R_NilValue,
                   int ncores = 1) {
  // Get number of alleles and markers.
  unsigned int m = MP.n_slices - 1;
  unsigned int p = MP.n_cols;
  unsigned int n = MP.n_rows;
  // Define covs as intercept.
  arma::mat covs = ones<mat>(n, 1);
  if (size_param.isNotNull()) {
    // Add other covariates.
    covs.insert_cols(1, Rcpp::as<arma::mat>(size_param));
  }
  // Get number of covariates (including intercept).
  unsigned int nCov = covs.n_cols;
  // Remove reference allele. Ref-1 because of 0-indexing.
  MP.shed_slice(ref - 1);
  // Compute inverse of t(M).
  arma::mat Mt = chol(sigma, "lower").i();
  // pre-multiply the IBD probablities with t(M).
  MP.each_slice([Mt](mat& X){ X = Mt * X; });
  // Pre-multiply the phenotype (y) with t(M).
  arma::mat tMy = Mt * y;
  // pre-multiply the intercept and covariates with t(M).
  arma::mat tMfixCovs = Mt * covs;
  // Compute residuals and RSS over all markers.
  arma::mat Vinv = inv_sympd(sigma);
  arma::mat Pr = Vinv - Vinv * covs * solve(covs.t() * Vinv * covs,
                                            covs.t() * Vinv);
  double RSSEnv = as_scalar(y.t() * Pr * y);
  // Define matrices for storing output.
  arma::mat beta1 = mat(tMfixCovs.n_cols, p);
  arma::mat beta2 = zeros<mat>(m, p);
  std::vector<double> FVal(p);
  std::vector<int> df1(p);
  std::vector<int> df2(p);
  // Loop over markers and compute betas and RSS.
#pragma omp parallel for num_threads(ncores)
  for (unsigned int i = 0; i < p; i ++) {
    arma::mat X = MP(span(), span(i), span());
    // Get indices of alleles in X that are not entirely 0.
    arma::uvec posInd = find( all(X) );
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
    arma::vec beta2j = tX2VinvX2Inv * tX2VinvY -
      tX2VinvX2Inv * tX1VinvX2.t() * beta1j;
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
    df1[i] = beta2j.n_elem;
    df2[i] = n - df1[i] - nCov;
    FVal[i] = (RSSEnv - RSSFullj) / RSSFullj * df2[i] / df1[i];
  }
  // Create output.
  List res;
  res["beta1"] = beta1;
  res["beta2"] = beta2;
  res["FVal"] = FVal;
  res["df1"] = df1;
  res["df2"] = df2;
  return res;
}
