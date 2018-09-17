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
List updatePrecCPP(unsigned int m,
                   unsigned int nc,
                   arma::mat omega,
                   arma::mat w,
                   arma::mat p,
                   bool het,
                   double maxDiag) {
  arma::mat pNew(size(p));
  arma::mat wNew(size(p));
  arma::mat cNew(size(p));
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
  return List::create(_["wNew"] = wNew,
                      _["pNew"] = pNew,
                      _["cNew"] = cNew);
}


