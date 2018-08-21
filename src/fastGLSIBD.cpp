#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List fastGLSIBDCPP(arma::cube MP,
                   arma::vec y,
                   arma::mat sigma,
                   arma::mat covs,
                   int ref) {
  // Get number of alleles, markers and covariates.
  int m = MP.n_slices - 1;
  int p = MP.n_cols;
  int n = MP.n_rows;
  int nCov = covs.n_cols;
  // Remove reference allele. Ref - 1 because of 0-indexing.
  MP.shed_slice(ref - 1);
  // Compute M.
  arma::mat M = (arma::chol(sigma)).i();
  // pre-multiply the IBD probablities with t(M).
  arma::cube tMPr = M.t() * MP.each_slice();
  // Pre-multiply the phenotype (y) with t(M).
  arma::mat tMy = M.t() * y;
  // pre-multiply the intercept and covariates with t(M).
  arma::mat tMfixCovs = M.t() * covs;
  // Compute residuals and RSS over all markers.
  mat Vinv = sigma.i();
  mat Pr = Vinv - Vinv * covs * solve(covs.t() * Vinv * covs, covs.t() * Vinv);
  double RSSEnv = as_scalar(y.t() * Pr * y);
  // Define matrices for storing output.
  arma::mat beta1 = mat(tMfixCovs.n_cols, p);
  arma::mat beta2 = zeros<mat>(m, p);
  std::vector<double> RSSFull(p);
  std::vector<double> FVal(p);
  std::vector<int> df1(p);
  std::vector<int> df2(p);
  for(int i = 0; i < p; i ++) {
    arma::mat X = tMPr( span(), span(i), span() );
    // Get indices of alleles in X that are not entirely 0.
    arma::uvec posInd = find(all(X) != 0 );
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
    RSSFull[i] = RSSFullj;
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







