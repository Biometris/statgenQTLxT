#include "EMFA.h"

#include <RcppArmadillo.h>
#include <testthat.h>

using namespace Rcpp;

context("EMFA Helper functions") {
  arma::mat c = "2 1; 1 2";
  arma::mat omega = "0.5 0.3; 0.3 0.5";
  arma::mat w = arma::mat(2, 1);
  arma::mat p = arma::mat(2, 2);
  updateFAHomVar(c, w, p, 1, 1e5);
  arma::mat wNew;
  arma::mat pNew;
  arma::mat cNew;
  updatePrec(1, 2, omega, w, p, wNew, pNew, cNew, false, 1e5);
  arma::mat wNew1;
  arma::mat pNew1;
  arma::mat cNew1;
  updatePrec(0, 2, omega, w, p, wNew1, pNew1, cNew1, false, 1e5);

  test_that("updatePrec produces correct output structure") {
    expect_true(wNew.n_rows == 2);
    expect_true(wNew.n_cols == 1);
    expect_true(pNew.n_rows == 2);
    expect_true(pNew.n_cols == 2);
    expect_true(cNew.n_rows == 2);
    expect_true(cNew.n_cols == 2);
    expect_true(wNew1.n_rows == 0);
    expect_true(approx_equal(cNew1, pNew1, "absdiff", 1e-8));
  }

  test_that("updatePrec functions properly") {
    arma::vec wRes = "0.547722557505166 0.547722557505166";
    arma::vec pRes = "5 0 0 5";
    arma::vec cRes = "3.125 -1.875 -1.875 3.125";
    expect_true(approx_equal(vectorise(wNew), wRes, "absdiff", 1e-8));
    expect_true(approx_equal(vectorise(pNew), pRes, "absdiff", 1e-8));
    expect_true(approx_equal(vectorise(cNew), cRes, "absdiff", 1e-8));
    arma::vec pRes1 = "2 0 0 2";
    expect_true(approx_equal(vectorise(pNew1), pRes1, "absdiff", 1e-8));
  }

  test_that("option het in updatePrec functions properly") {
    arma::mat wNewa;
    arma::mat pNewa;
    arma::mat cNewa;
    updatePrec(1, 2, omega, w, p, wNewa, pNewa, cNewa, true, 1e5);
    arma::vec wResa = "-2.53587580354153e-05 2.53587580354153e-05";
    arma::vec pResa = "2.00000000257227 0 0 2.00000000257227";
    arma::vec cResa = "2 2.57226643639495e-09 2.57226643639495e-09 2";
    arma::mat wNew1a;
    arma::mat pNew1a;
    arma::mat cNew1a;
    updatePrec(0, 2, omega, w, p, wNew1a, pNew1a, cNew1a, true, 1e5);
    arma::vec pRes1a = "2 0 0 2";
    expect_true(approx_equal(vectorise(pNew1a), pRes1a, "absdiff", 1e-8));
  }

  test_that("vecInvDiag functions properly") {
    arma::mat res = "4 7; 5 9; 6 11;";
    res = 1 / res;
    arma::mat vInvDiag = vecInvDiag(arma::vec {1, 2}, arma::vec {3, 4, 5});
    expect_true(vInvDiag.n_cols == 2);
    expect_true(vInvDiag.n_rows == 3);
    expect_true(approx_equal(vInvDiag, res, "absdiff", 1e-8));
  }

  test_that("tracePInvDiag functions properly") {
    arma::vec res = "0.616666666666667, 0.344877344877345";
    arma::vec tpid = tracePInvDiag(arma::vec {1, 2}, arma::vec {3, 4, 5});
    expect_true(tpid.n_elem == 2);
    expect_true(approx_equal(tpid, res, "absdiff", 1e-8));
  }
}
