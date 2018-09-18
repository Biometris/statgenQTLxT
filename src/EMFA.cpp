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
void updateFAHomVarCPP(arma::mat s,
                       arma::mat& wNew,
                       arma::mat& pNew,
                       unsigned int m,
                       double maxDiag = 1e4) {
  unsigned int nc = s.n_cols;
  arma::vec eigVals(nc);
  arma::mat eigVecs(nc, nc);
  eig_sym(eigVals, eigVecs, s);
  eigVals = reverse(eigVals);
  eigVecs = fliplr(eigVecs);
  double sigma2 = mean(eigVals.tail(nc - m));
  sigma2 = std::max(sigma2, 1 / maxDiag);
  // Split cases for extra robustness.
  if (m == 1) {
    wNew = eigVecs.head_cols(1) * sqrt(eigVals.head(1) - sigma2);
  } else {
    wNew = eigVecs.head_cols(m) * diagmat(sqrt(eigVals.head(m) - sigma2));
  }
  pNew = eye<mat>(nc, nc) / sigma2;
}


// [[Rcpp::export]]
void updateFACPP(arma::mat y,
                 arma::mat wStart,
                 arma::mat pStart,
                 arma::mat& wNew,
                 arma::mat& pNew,
                 unsigned int m0,
                 bool hetVar = false,
                 double maxDiag = 1e4,
                 double tolerance = 1e-4,
                 unsigned int maxIter = 100,
                 bool printProgress = false) {
  unsigned int nc = y.n_cols;
  unsigned int nr = y.n_rows;
  // Set start values for P and W.
  arma::mat w(nc, nc);
  arma::mat p(nc, nc);
  unsigned int m = 2;
  if (m0) {
    m = m0;
  }
  if (wStart.n_elem > 0) {
    w = wStart;
    p = pStart;
    if (!m0) {
      m = w.n_cols;
    }
  } else {
    arma::vec eigVals(nc);
    arma::mat eigVecs(nc, nc);
    eig_sym(eigVals, eigVecs, y * y.t() / nr);
    eigVals = reverse(eigVals);
    eigVecs = fliplr(eigVecs);
    double sigma2 = mean(eigVals.tail(nc - m));
    p = eye<mat>(nc, nc) / sigma2;
    w = eigVecs.head_cols(m) * diagmat(sqrt(eigVals.head(m) - sigma2));
  }
  // Set start values for iterations and difference.
  double totDiff = tolerance + 1;
  // EM
  wNew = mat(size(w));
  pNew = p;
  for (unsigned int i = 0; i < maxIter && totDiff > tolerance; i++) {
    if (m == 1) {
      arma::mat m1(m, nc);
      arma::mat sigma(m, m);
      if (hetVar) {
        arma::mat b = w.t() * p * w;
        sigma = 1 / (1 + b);
        m1 = as_scalar(sigma) * w.t() * p * y.t();
      } else {
        arma::mat b = p(0, 0) * w.t() * w;
        sigma = 1 / (1 + b);
        m1 = p(0, 0) * as_scalar(sigma) * (y * w).t();
      }
      double a = 1 / (nr * as_scalar(sigma) + accu(square(m1)));
      wNew = a * y.t() * m1.t();
      if (hetVar) {
        arma::rowvec d1 = sum(square(y));
        arma::vec d2 = square(vectorise(wNew)) / a;
        arma::mat d3 = wNew * m1 * y;
        arma::vec dTot = d1.t() + d2 - 2 * d3.diag();
        pNew.diag() = nr / dTot;
      } else {
        arma::mat dNew = y.t() - wNew * m1;
        arma::mat sNew = dNew * dNew.t() / nr;
        pNew.diag().fill( 1 / mean(sNew.diag()) );
      }
    } else {
      arma::mat m1(m, nc);
      arma::mat sigma(m, m);
      if (hetVar) {
        arma::mat b = w.t() * p * w;
        sigma = (b + eye<mat>(m, m)).i();
        m1 = sigma * w.t() * p * y.t();
      } else {
        arma::mat b = p(0, 0) * w.t() * w;
        sigma = (b + eye<mat>(m, m)).i();
        m1 = p(0, 0) * sigma * (y * w).t();
      }
      arma::mat a = (nr * sigma + m1 * m1.t());
      wNew = y.t() * m1.t() * a.i();
      if (hetVar) {
        arma::rowvec d1 = sum(square(y));
        arma::mat d2 = (wNew * a * wNew.t());
        arma::mat d3 = wNew * m1 * y;
        arma::vec dTot = d1.t() + d2.diag() - 2 * d3.diag();
        pNew.diag() = nr / dTot;
      } else {
        arma::mat dNew = y.t() - wNew * m1;
        arma::mat sNew = dNew * dNew.t() / nr;
        pNew.diag().fill( 1 / mean(sNew.diag()) );
      }
    }
    pNew.diag() = clamp(diagvec(pNew), min(diagvec(pNew)), maxDiag);
    double pDiff = accu(abs(pNew - p));
    double wDiff = accu(abs(wNew - w));
    totDiff = pDiff + wDiff;
    if (printProgress) {
      Rcout << "Iteration " + std::to_string(i + 1) + " : " +
        std::to_string(pDiff) + "  " + std::to_string(wDiff) << endl;
    }
    // Set values for next iteration
    p = pNew;
    w = wNew;
  }
}

// [[Rcpp::export]]
void updatePrecCPP(unsigned int m,
                   unsigned int nc,
                   arma::mat omega,
                   const arma::mat w,
                   const arma::mat p,
                   arma::mat& wNew,
                   arma::mat& pNew,
                   arma::mat& cNew,
                   bool het,
                   double maxDiag) {
  if (m == 0) {
    // Recall that the model is C^{-1} = P^{-1} + W W^t.
    // when m == 0, W = 0 and Cm = P.
    if (nc > 1) {
      if (het) {
        arma::vec tau = 1 / omega.diag();
        pNew = diagmat(clamp(tau, tau.min(), maxDiag));
      } else {
        double tau = std::min(maxDiag, nc / accu(omega.diag()));
        pNew <- tau * eye<mat>(nc, nc);
      }
    } else {
      pNew <- 1 / omega;
    }
    cNew = pNew;
  } else {
    // When rank(Omega) = Q, A should be the Q x p matrix such that
    // Omega = A^t A / Q.
    arma::mat a = sqrtmat_sympd(omega) * sqrt(omega.n_rows);
    if (het) {
      updateFACPP(a, w, p, wNew, pNew, w.n_cols, het, maxDiag);
    } else {
      updateFAHomVarCPP(omega, wNew, pNew, m);
    }
    cNew = (pNew.i() + wNew * wNew.t()).i();
  }
}

arma::mat vecInvDiagCPP(arma::vec x,
                        arma::vec y) {
  arma::mat z = zeros<mat>(y.n_elem, x.n_elem);
  for (unsigned int i = 0; i < z.n_cols; i++) {
    z.col(i) = 1 / (1 + x(i) * y);
  }
  return z;
}

arma::vec tracePInvDiagCPP(arma::vec x,
                           arma::vec y) {
  arma::vec z = zeros<vec>(x.n_elem);
  for (unsigned int i = 0; i < z.n_elem; i++) {
    z(i) = accu(1 / (1 + x(i) * y));
  }
  return z;
}

// [[Rcpp::export]]
List EMFACPP(arma::mat y,
             arma::mat k,
             Rcpp::Nullable<Rcpp::NumericVector> size_param_x = R_NilValue,
             bool cmHet = false,
             bool dmHet = false,
             double tolerance = 1e-4,
             unsigned int maxIter = 300L,
             Rcpp::Nullable<Rcpp::NumericVector> size_param_cmStart = R_NilValue,
             Rcpp::Nullable<Rcpp::NumericVector> size_param_dmStart = R_NilValue,
             unsigned int mG = 1,
             unsigned int mE = 1,
             double maxDiag = 1e4,
             bool stopIfDecreasing = false) {
  arma::mat x;
  if (size_param_x.isNotNull()) {
    x = Rcpp::as<arma::mat>(size_param_x);
  } else {
    x = ones<mat>(k.n_rows, 1);
  }
  unsigned int nc = x.n_cols;
  unsigned int n = k.n_cols;
  unsigned int p = y.n_cols;
  arma::mat b = zeros<mat>(nc, p);
  arma::mat xtxInvxt = solve(x.t() * x, x.t());
  arma::vec dk(n);
  arma::mat uk(n, n);
  eig_sym(dk, uk, inv_sympd(k));
  dk = reverse(dk);
  uk = fliplr(uk);
  arma::mat lambdaR = diagmat(dk);
  // Set starting values for Cm.
  arma::mat cm;
  if (size_param_cmStart.isNotNull()) {
    cm = Rcpp::as<arma::mat>(size_param_cmStart);
  } else {
    if (mG == 0) {
      cm = eye<mat>(p, p) * 2;
    } else {
      cm = inv_sympd(cor(y) + eye<mat>(p, p) / 4);
    }
  }
  // Set starting values for Dm.
  arma::mat dm;
  if (size_param_dmStart.isNotNull()) {
    cm = Rcpp::as<arma::mat>(size_param_dmStart);
  } else {
    if (mE == 0) {
      dm = eye<mat>(p, p) * 2;
    } else {
      dm = inv_sympd(cor(y) + eye<mat>(p, p) / 4);
    }
  }
  // The model is Cm^{-1} = P^{-1} + W W^t
  // Given a starting value for Cm, set starting values for P and W
  arma::mat wg(nc, p);
  arma::mat pg(p, p);
  updateFAHomVarCPP(inv_sympd(cm), wg, pg, mG);
  // Given a starting value for Dm, set starting values for P and W
  arma::mat we(nc, p);
  arma::mat pe(p, p);
  updateFAHomVarCPP(inv_sympd(dm), we, pe, mE);
  // Set further starting values.
  double eLogLik = 0;
  arma::mat mu = zeros<mat>(n, p);
  // EM following the notation of Dahl et al.
  for (unsigned int i = 0; i < maxIter; i++) {
    arma::mat dmSqrt = sqrtmat_sympd(dm);
    arma::mat dmSqrtInv = inv_sympd(dmSqrt);
    arma::vec v1(p);
    arma::mat q1(p, p);
    eig_sym(v1, q1, dmSqrtInv * cm * dmSqrtInv);
    v1 = reverse(v1);
    q1 = fliplr(q1);
    arma::mat cmSqrt = sqrtmat_sympd(cm);
    arma::mat cmSqrtInv = inv_sympd(cmSqrt);
    arma::vec v2(p);
    arma::mat q2(p, p);
    eig_sym(v2, q2, cmSqrtInv * dm * cmSqrtInv);
    v2 = reverse(v2);
    q2 = fliplr(q2);
    arma::mat s1(size(y));
    arma::mat s2(size(y));
    if (nc > 0) {
      arma::mat tUyminxb = uk.t() * (y - x * b);
      s1 = vecInvDiagCPP(v1, dk) % (tUyminxb * dmSqrt * q1);
      s2 = vecInvDiagCPP(v2, 1 / dk) % (tUyminxb * cmSqrt * q2);
    } else {
      s1 = vecInvDiagCPP(v1, dk) % (uk.t() * y * dmSqrt * q1);
      s2 = vecInvDiagCPP(v2, 1 / dk) % (uk.t() * y * cmSqrt * q2);
    }
    arma::vec trP1 = tracePInvDiagCPP(v1, dk);
    arma::vec trP2 = tracePInvDiagCPP(v2, 1 / dk);
    arma::mat part1 = zeros<mat>(p, p);
    arma::mat part2 = zeros<mat>(p, p);
    if (p > 1) {
      part1 = dmSqrtInv * q1 * diagmat(trP1) * q1.t() * dmSqrtInv;
      part2 = cmSqrtInv * q2 * diagmat(trP2) * q2.t() * cmSqrtInv;
    } else {
      part1 = dmSqrtInv * q1 * as_scalar(trP1) * q1.t() * dmSqrtInv;
      part2 = cmSqrtInv * q2 * as_scalar(trP2) * q2.t() * cmSqrtInv;
    }
    arma::mat part3 = (cmSqrtInv * q2 * s2.t()) * (cmSqrtInv * q2 * s2.t()).t();
    arma::mat part4 = dmSqrtInv * q1 * s1.t() * lambdaR *
      s1 * (dmSqrtInv * q1).t();
    if (nc > 0) {
      mu = uk * s1 * (dmSqrtInv * q1).t();
      b = xtxInvxt * (y - mu);
    }
    arma::mat omega1 = (part1 + part3) / n;
    arma::mat omega2 = (part2 + part4) / n;
    // Update Cm
    updatePrecCPP(mG, p, omega2, wg, pg, wg, pg, cm, cmHet, maxDiag);
    // Update Dm
    updatePrecCPP(mE, p, omega1, we, pe, we, pe, dm, dmHet, maxDiag);
    // Compute log-likelihood and check stopping criteria.
    double eLogLikOld = eLogLik;
    double eLogLikCm = n * (log_det(cm).real() - trace(cm * omega2));
    double eLogLikDm = n * (log_det(dm).real() - trace(dm * omega1));
    eLogLik = eLogLikCm + eLogLikDm;
    if (stopIfDecreasing && i > 49 && eLogLik < eLogLikOld - 0.1) {
      break;
    }
    if (i > 0 && (fabs(eLogLik - eLogLikOld) < tolerance)) {
      break;
    }
  }
  return List::create(_["vg"] = inv_sympd(cm),
                      _["ve"] = inv_sympd(dm));
}




