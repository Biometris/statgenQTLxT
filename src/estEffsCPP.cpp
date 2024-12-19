#include <RcppArmadillo.h>
#include "getThr.h"
#include "utils.h"

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
//' @inheritParams estEffsCPP
//'
//' @param size_param_x An optional c x n covariate matrix, c being the number
//' of covariates and n being the number of genotypes. c has to be at least one
//' (typically an intercept). No missing values are allowed.
//' @param vInv A n x p x p cube containing for each genotype l the
//' p x p matrix \eqn{v_l ^ {-1}} (in the notation of Zhou and Stephens).
//'
//' @returns A numerical value for the \eqn{t(y) * P * y} part of the
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
                         Rcpp::Nullable<Rcpp::NumericVector> size_param_x = R_NilValue) {
  unsigned int nc = 0;
  arma::mat x;
  // Define covs as intercept.
  if (size_param_x.isNotNull()) {
    // Add other covariates.
    x = Rcpp::as<arma::mat>(size_param_x);
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
//' @param y0 An n x p matrix of observed phenotypes, on p traits or environments
//' for n genotypes. No missing values are allowed.
//' @param w0 An n x c covariate matrix, c being the number of covariates and n
//' being the number of genotypes. c has to be at least one (typically an
//' intercept). No missing values are allowed.
//' @param x0 An n x ns matrix of marker scores. Neither missing values nor
//' non-segregating markers are allowed.
//' @param vg A p x p matrix of genetic covariances.
//' @param ve A p x p matrix of environmental covariances.
//' @param k An n x n genetic relatedness matrix.
//' @param returnSe Should standard errors and p-values be returned?
//' @param estCom Should the common SNP-effect model be fitted?
//' @param nCores An integer indicating the number of cores used for parallel
//' computation.
//'
//' @returns A list containing the estimates, optionally the standard errors of
//' the estimates and corresponding p-values. If \code{estCom = TRUE} also
//' common SNP-effects, their standard errors and corresponding p-values and
//' the p-values for QtlxE are output.
//'
//' @references Zhou, X. and Stephens, M. (2014). Efficient multivariate linear
//' mixed model algorithms for genome-wide association studies. Nature Methods,
//' February 2014, Vol. 11, p. 407–409
//' @keywords internal
// [[Rcpp::export]]
List estEffsCPP(arma::mat& y0,
                arma::mat& w0,
                arma::mat& x0,
                arma::mat& vg,
                arma::mat& ve,
                arma::mat& k,
                bool returnSe = true,
                bool estCom = false,
                Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  // Extract number of traits, individuals, covariates and SNPs.
  unsigned int p = y0.n_cols;
  unsigned int n = y0.n_rows;
  unsigned int nc = w0.n_cols;
  unsigned int ns = x0.n_cols;
  unsigned int nChunks = R::fmax2(R::fmin2(ns, ceil(p * nc * ns / 50000)), 1);
  k = nearestPDQTLxT(k, false, false, true, false, true, 1e-10, 1e-7, 1e-8, 100);
  arma::mat uk;
  arma::vec dk;
  arma::eig_sym(dk, uk, k);
  // Transform y, w and x.
  arma::mat y = y0.t() * uk;
  arma::mat w = w0.t() * uk;
  arma::mat x = x0.t() * uk;
  //Square each element of X.
  arma::mat x2 = square(x);
  arma::cube vInv = cube(p, p, n);
  arma::mat vInvy = mat(p, n);
  int nThr = getThr(nCores);
#pragma omp parallel for num_threads(nThr)
  for(unsigned int i = 0; i < n; i++) {
    //vInv.slice(i) = inv_sympd(dk(i) * vg + ve);
    vInv.slice(i) = (dk(i) * vg + ve).i();
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
    for (unsigned int cv = 0; cv < nc; cv++) {
      s3.each_row() %= w.row(cv);
    }
  }
  arma::mat eff = zeros<mat>(p, ns);
  arma::mat effSe = zeros<mat>(p, ns);
  NumericVector pVals = NumericVector(ns);
  NumericVector effCom = NumericVector(ns);
  NumericVector effSeCom = NumericVector(ns);
  NumericVector pValsCom = NumericVector(ns);
  NumericVector pValsQtlE = NumericVector(ns);
  unsigned int nsCh = ns / nChunks;
  // Loop over chunks.
  for (unsigned int ch = 0; ch < nChunks; ch++) {
    int chStart = ch * nsCh;
    if (ch == nChunks - 1 && ns % nsCh > 0) {
      nsCh += ns % nsCh;
    }
    int chEnd = chStart + nsCh - 1;
    // Define SNP-dependent quantities.
    arma::cube x2Vinvx1 = zeros<cube>(nc, nsCh, p * p);
    // Fill X2VinvX1Arr by looping over covariates.
    for (unsigned int cv = 0; cv < nc; cv++) {
      arma::mat vInvMatCv = vInvMat;
      vInvMatCv.each_row() %= w.row(cv);
      x2Vinvx1( span(cv), span(), span() ) = x.rows(chStart, chEnd) * vInvMatCv.t();
    }
    // Compute VBeta and V per SNP - computation only depends on individuals.
    arma::mat vBetaSnp = x2.rows(chStart, chEnd) * vInvMat.t();
    arma::mat vSnp = x.rows(chStart, chEnd) * vInvy.t();
    arma::mat vSnpQ = zeros(nsCh, p);
    if (returnSe) {
      // Compute VSnpQ by looping over individuals.
      for (unsigned int i = 0; i < n; i++) {
        vSnpQ += x.rows(chStart, chEnd).col(i) * vInvy.col(i).t();
      }
    }
    arma::vec ss1 = zeros<vec>(nsCh);
    // Compute degrees of freedom for full model.
    double dfFull = (n - nc - 1) * p;
    // The remaining calculations are SNP-dependent.
#pragma omp parallel for num_threads(nThr)
    for (unsigned int snp = 0; snp < nsCh; snp++) {
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
      eff.col(snp + chStart) = vBetaSnpInv *
        (vSnp.row(snp).t() - x2Vinvx1Snp.t() * effCovSnp);
      if (returnSe) {
        // Compute inverse of QSnp.
        arma::mat qSnpInv = (join_cols(join_rows(vBeta, x2Vinvx1Snp),
                                       join_rows(x2Vinvx1Snp.t(),
                                                 vBetaSnpInv.i()))).i();
        // Compute SE of SNP effect.
        arma::vec effSeSnp = sqrt(qSnpInv.diag());
        effSe.col(snp + chStart) = effSeSnp.tail(p);
        // Compute SS1 per SNP.
        arma::mat qSnp = join_cols(vq, vSnpQ.row(snp).t());
        ss1(snp) = qScal - as_scalar(qSnp.t() * qSnpInv * qSnp);
        // Compute p-values for SNP effect.
        double fValx = ((ss0 - ss1(snp)) / ss1(snp)) * dfFull / p;
        pVals(snp + chStart) = R::pf(fValx, p, dfFull, false, false);
      } // End returnSe
    } // End loop over SNPs
    arma::vec ss1Com = zeros<vec>(nsCh);
    // Compute degrees of freedom for common effect model.
    double dfCom  = (n - nc) * p - 1;
    if (estCom) {
      // Calculations for common effect are similar to those above.
      // However the dimensions are lower and therefore some calculations
      // have been simplified.
      // Define SNP-dependent quantities.
      arma::cube x2Vinvx1Com = cube(nc, nsCh, p);
      // Fill x2Vinvx1Com by looping over covariates.
      for (unsigned int cv = 0; cv < nc; cv++) {
        x2Vinvx1Com( span(cv), span(), span() ) = x.rows(chStart, chEnd) * s3.t();
      }
      // Compute common VBeta and V per SNP - only depends on individuals.
      arma::mat vSnpCom = x.rows(chStart, chEnd) * s1.t();
      arma::vec vBetaSnpCom = x2.rows(chStart, chEnd) * s2.t();
      arma::colvec vSnpQCom;
      if (returnSe) {
        // Compute VSnpQ for common effect.
        vSnpQCom = sum(vSnpCom, 1);
      }
      // The remaining calculations are SNP-dependent.
#pragma omp parallel for num_threads(nThr)
      for (unsigned int snp = 0; snp < nsCh; snp++) {
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
        effCom(snp + chStart) = as_scalar(vBetaSnpInvCom *
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
          effSeCom(snp + chStart) = as_scalar(effSeComSnp.tail(1));
          // Compute SS1 for common effect.
          arma::vec qSnpCom = vq;
          qSnpCom.resize(qSnpCom.size() + 1);
          qSnpCom(qSnpCom.size() - 1) = vSnpQCom(snp);
          ss1Com(snp) = qScal - as_scalar(qSnpCom.t() * qSnpInvCom * qSnpCom);
          // Compute p-values for common SNP effect.
          double fValComx = (ss0 - ss1Com(snp)) / ss1Com(snp) * dfCom;
          pValsCom(snp + chStart) = R::pf(fValComx, 1, dfCom, false, false);
          double fValQtlEx = ((ss1Com(snp) - ss1(snp)) / ss1(snp)) *
            dfFull / (p - 1);
          pValsQtlE(snp + chStart) = R::pf(fValQtlEx, p - 1, dfFull, false, false);
        } // End returnSe
      } // End loop over SNPs
    } // End estCom
  } // End loop over chuncks
  // Construct output.
  return List::create(_["effs"] = eff,
                      _["effsSe"] = effSe,
                      _["pVals"] = pVals,
                      _["effsCom"] = effCom,
                      _["effsComSe"] = effSeCom,
                      _["pValsCom"] = pValsCom,
                      _["pValsQtlE"] = pValsQtlE);
}
