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

//' Compute tYPY as in Zhou and Stephens eqn. 50.
//'
//' Compute \eqn{t(y) * P * y}, part of the log-likelihood functions from
//' equation 26 and 27 in Zhou and Stephens using equation 50. Equation 56, 57
//' and 58 are used to do the actual computations.
//'
//' It is assumed that X and Y have already been rotated by Uk, where Uk is such
//' that the kinship matrix K equals \eqn{K = Uk * Dk * t(Uk)}.\cr
//' The original X and Y are right multiplied by Uk, e.g. \code{Y <- Y * Uk}.
//' See Zhou and Stephens (2014), supplement.\cr
//' It is these rotated versions that are the input of this function.
//'
//' @inheritParams LLDiag
//'
//' @return A numerical value for the \eqn{t(y) * P * y} part of the
//' log-likelihood function.
//'
//' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
//' mixed model algorithms for genome-wide association studies. Nature Methods,
//' February 2014, Vol. 11, p. 407–409
//'
//' @keywords internal
// [[Rcpp::export]]
double LLQuadFormDiagCPP(const arma::mat &y,
                         const arma::cube &vInv,
                         Rcpp::Nullable<Rcpp::NumericVector> size_param = R_NilValue) {
  unsigned int nc = 0;
  arma::mat x;
  // Define covs as intercept.
  if (size_param.isNotNull()) {
    // Add other covariates.
    x = Rcpp::as<arma::mat>(size_param);
    nc = x.n_rows;
  }
  // Get number of genotypes and markers.
  unsigned int n = y.n_cols;
  unsigned int p = y.n_rows;
  // Compute scalair part.
  double qScal = 0;
  for (unsigned int i = 0; i < n; i++) {
    qScal += as_scalar(y.col(i).t() * vInv.slice(i) * y.col(i));
  }
  if (nc > 0) {
    arma::vec qVec = zeros<vec>(p);
    arma::mat qMat = zeros<mat>(p, p * nc);
    // compute q and Q.
    for (unsigned int i = 0; i < n; i++) {
      qVec += kron(x.col(i), vInv.slice(i) * y.col(i));
      qMat += kron(x.col(i) * x.col(i).t(), vInv.slice(i));
    }
    // Compute quadratic part.
    arma::mat qQuad = qVec.t() * solve(qMat, qVec);
    qScal -= as_scalar(qQuad * qQuad.t());
  }
  return qScal;
}

//' Estimates for covariates
//'
//' Compute the estimates and standard errors for the covariates in the input
//' matrix W.
//'
//' @param y An n x p matrix of observed phenotypes, on p traits or environments
//' for n genotypes. No missing values are allowed.
//' @param w An n x c covariate matrix, c being the number of covariates and n
//' being the number of genotypes. c has to be at least one (typically an
//' intercept). No missing values are allowed.
//' @param x An n x ns matrix of marker scores. Neither missing values nor
//' non-segregating markers are allowed.
//' @param vg A p x p matrix of genetic covariances.
//' @param ve A p x p matrix of environmental covariances.
//' @param k An n x n genetic relatedness matrix.
//' @param returnSe Should standard errors and p-values be returned?
//' @param estCom Should the common SNP-effect model be fitted?
//' @param ncores An integer indicating the number of cores used for parallel
//' computation.
//'
//' @return A list containing the estimates, optionally the standard errors of
//' the estimates and corresponding p-values. If \code{estCom = TRUE} also
//' common SNP-effects, their standard errors and corresponding p-values and
//' the p-values for QtlxE are output.
//'
//' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
//' mixed model algorithms for genome-wide association studies. Nature Methods,
//' February 2014, Vol. 11, p. 407–409
// [[Rcpp::export]]
List estEffsCPP(arma::mat y,
                arma::mat w,
                arma::mat x,
                arma::mat& vg,
                arma::mat& ve,
                arma::mat& k,
                bool returnSe = true,
                bool estCom = false,
                int ncores = 4) {
  // Extract number of traits, individuals, covariates and SNPs.
  unsigned int p = y.n_cols;
  unsigned int n = y.n_rows;
  unsigned int nc = w.n_cols;
  unsigned int ns = x.n_cols;
  arma::mat uk;
  arma::vec dk;
  arma::eig_sym(dk, uk, k);
  dk = clamp(dk, dk.max() * pow(10, -5), dk.max());
  // Transform y, w and x.
  y = y.t() * uk;
  w = w.t() * uk;
  x = x.t() * uk;
  //Square each element of X.
  arma::mat x2 = square(x);
  arma::cube vInv = cube(p, p, n);
  arma::mat vInvy = mat(p, n);
#pragma omp parallel for num_threads(ncores)
  for(unsigned int i = 0; i < n; i++) {
    vInv.slice(i) = inv_sympd(dk(i) * vg + ve);
    // vInvy is used is several equations so is computed once here.
    vInvy.col(i) = vInv.slice(i) * y.col(i);
  }
  // create a second instance of vInvLst, in matrix form.
  arma::mat vInvMat = reshape( mat(vInv.memptr(), vInv.n_elem, 1, false),
                               p * p, n );
  // Compute quantities that are independent of the SNPs.
  arma::vec vVec = zeros<vec>(p * nc);
  arma::mat vBeta = zeros<mat>(p * nc, p * nc);
  for (unsigned int i = 0; i < n; i++) {
    vVec += kron(w.col(i), vInvy.col(i));
    vBeta += kron(w.col(i) * w.col(i).t(), vInv.slice(i));
  }
  arma::vec vq;
  double qScal = 0;
  double ss0 = 0;
  if (returnSe) {
    // Compute SS0 for null model with the trait specific means only.
    arma::mat est0 = solve(vBeta, vVec);
    est0.reshape( p, est0.n_rows / p );
    arma::mat fitMean0 = est0 * w;
    ss0 = LLQuadFormDiagCPP(y - fitMean0, vInv);
    // Compute VQ
    vq = vectorise(vInvy * w.t());
    // Compute scalar part of SS1.
    qScal = accu(y % vInvy);
  }
  arma::rowvec s1, s2;
  arma::mat s3;
  if (estCom) {
    // Compute chunk independent quantities.
    s1 = sum(vInvy, 0);
    s3 = sum(vInv, 1);
    s2 = sum(s3, 0);
  }
  // Define SNP-dependent quantities.
  arma::cube x2Vinvx1 = zeros<cube>(nc, ns, p * p);
  // Fill X2VinvX1Arr by looping over covariates.
  for (unsigned int cv = 0; cv < nc; cv++) {
    arma::mat vInvMatCv = vInvMat;
    vInvMatCv.each_row() %= w.row(cv);
    x2Vinvx1( span(cv), span(), span() ) = x * vInvMatCv.t();
  }
  // Compute VBeta and V per SNP - computation only depends on individuals.
  arma::mat vBetaSnp = x2 * vInvMat.t();
  arma::mat vSnp = x * vInvy.t();
  arma::mat vSnpQ = zeros(ns, p);
  if (returnSe) {
    // Compute VSnpQ by looping over individuals.
    for (unsigned int i = 0; i < n; i++) {
      vSnpQ += x.col(i) * vInvy.col(i).t();
    }
  }
  arma::mat eff = zeros<mat>(p, ns);
  arma::mat effSe = zeros<mat>(p, ns);
  arma::vec ss1 = zeros<vec>(ns);
  arma::vec fVals = zeros<vec>(ns);
  // Compute degrees of freedom for full model.
  double dfFull = (n - nc - 1) * p;
  // The remaining calculations are SNP-dependent.
#pragma omp parallel for num_threads(ncores)
  for (unsigned int snp = 0; snp < ns; snp++) {
    // Compute inverse of VBetaSnp.
    arma::mat vBetaSnpInv = vBetaSnp(snp, span());
    vBetaSnpInv.reshape( p, p );
    vBetaSnpInv = vBetaSnpInv.i();
    // Extract X2VinvX1 from array.
    arma::mat x2Vinvx1Snp = x2Vinvx1( span(), span(snp), span() );
    x2Vinvx1Snp = x2Vinvx1Snp.t();
    x2Vinvx1Snp.reshape( p, p * nc);
    x2Vinvx1Snp = x2Vinvx1Snp.t();
    // Compute inverse of XS.
    arma::mat xsInv = (vBeta - x2Vinvx1Snp * vBetaSnpInv * x2Vinvx1Snp.t()).i();
    // Compute SNP dependent effect.
    arma::mat effCovSnp = xsInv *
      (vVec - x2Vinvx1Snp * vBetaSnpInv * vSnp.row(snp).t());
    // // Compute SNP effect.
    eff.col(snp) = vBetaSnpInv *
      (vSnp.row(snp).t() - x2Vinvx1Snp.t() * effCovSnp);
    if (returnSe) {
      // Compute inverse of QSnp.
      arma::mat qSnpInv = (join_cols(join_rows(vBeta, x2Vinvx1Snp),
                                     join_rows(x2Vinvx1Snp.t(),
                                               vBetaSnpInv.i()))).i();
      // Compute SE of SNP effect.
      arma::vec effSeSnp = sqrt(qSnpInv.diag());
      effSe.col(snp) = effSeSnp.tail(p);
      // Compute SS1 per SNP.
      arma::mat qSnp = join_cols(vq, vSnpQ.row(snp).t());
      ss1(snp) = qScal - as_scalar(qSnp.t() * qSnpInv * qSnp);
      // Compute F-values for SNP effect.
      fVals(snp) = ((ss0 - ss1(snp)) / ss1(snp)) * dfFull / p;
    } // End returnSe
  } // End loop over SNPs
  arma::vec effCom = zeros<vec>(ns);
  arma::vec effSeCom = zeros<vec>(ns);
  arma::vec ss1Com = zeros<vec>(ns);
  arma::vec fValsCom = zeros<vec>(ns);
  arma::vec fValsQtlE = zeros<vec>(ns);
  // Compute degrees of freedom for common effect model.
  double dfCom  = (n - nc) * p - 1;
  if (estCom) {
    // Calculations for common effect are similar to those above.
    // However the dimensions are lower and therefore some calculations
    // have been simplified.
    // Define SNP-dependent quantities.
    arma::cube x2Vinvx1Com = cube(nc, ns, p);
    // Fill X2VinvX1ArrCom by looping over covariates.
    for (unsigned int cv = 0; cv < nc; cv++) {
      s3.each_row() %= w.row(cv);
      x2Vinvx1Com( span(cv), span(), span() ) = x * s3.t();
    }
    // Compute common VBeta and V per SNP - only depends on individuals.
    arma::mat vSnpCom = x * s1.t();
    arma::vec vBetaSnpCom = x2 * s2.t();
    arma::colvec vSnpQCom;
    if (returnSe) {
      // Compute VSnpQ for common effect.
      vSnpQCom = sum(vSnpCom, 1);
    }
    // The remaining calculations are SNP-dependent.
#pragma omp parallel for num_threads(ncores)
    for (unsigned int snp = 0; snp < ns; snp++) {
      // Compute VBetaSnpInvCom - scalar value
      double vBetaSnpInvCom = 1 / vBetaSnpCom(snp);
      // Extract X2VinvX1Com from array.
      arma::vec x2Vinvx1ComSnp =
        vectorise(x2Vinvx1Com( span(), span(snp), span() ));
      // Compute inverse of XS for common effect.
      arma::mat xsInvCom = (vBeta -
        vBetaSnpInvCom * x2Vinvx1ComSnp * x2Vinvx1ComSnp.t()).i();
      // Compute SNP dependent common effect.
      arma::mat effCovComSnp = xsInvCom *
        (vVec - vBetaSnpInvCom * vSnpCom(snp) * x2Vinvx1ComSnp);
      // // Compute common effect.
      effCom(snp) = as_scalar(vBetaSnpInvCom *
        (vSnpCom(snp) - x2Vinvx1ComSnp.t() * effCovComSnp));
      if (returnSe) {
        // Compute inverse of QSnp for common effect.
        arma::vec lastRow = x2Vinvx1ComSnp;
        lastRow.resize(lastRow.size() + 1);
        lastRow(lastRow.size() - 1) = 1 / vBetaSnpInvCom;
        arma::mat qSnpInvCom = (join_cols(join_rows(vBeta, x2Vinvx1ComSnp),
                                          lastRow.t())).i();
        // Compute SE of common SNP effect.
        arma::vec effSeComSnp = sqrt(qSnpInvCom.diag());
        effSeCom(snp) = as_scalar(effSeComSnp.tail(1));
        // Compute SS1 for common effect.
        arma::vec qSnpCom = vq;
        qSnpCom.resize(qSnpCom.size() + 1);
        qSnpCom(qSnpCom.size() - 1) = vSnpQCom(snp);
        ss1Com(snp) = qScal - as_scalar(qSnpCom.t() * qSnpInvCom * qSnpCom);
        // Compute F-values for common SNP effect.
        fValsCom(snp) = (ss0 -  ss1Com(snp)) /  ss1Com(snp) * dfCom;
        fValsQtlE(snp) = ((ss1Com(snp) - ss1(snp)) / ss1(snp)) * dfFull / (p - 1);
      } // End returnSe
    } // End loop over SNPs
  } // End estCom
  // Construct output.
  return List::create(_["effs"] = eff,
                      _["effsSe"] = effSe,
                      _["fVals"] = fVals,
                      _["effsCom"] = effCom,
                      _["effsComSe"] = effSeCom,
                      _["fValCom"] = fValsCom,
                      _["fValQtlE"] = fValsQtlE,
                      _["dfFull"] = dfFull,
                      _["dfCom"] = dfCom);

}
