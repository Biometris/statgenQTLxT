#ifndef EMFA_H_
#define EMFA_H_

#include "RcppArmadillo.h"

Rcpp::List EMFA(arma::mat y,
                arma::mat k,
                Rcpp::Nullable<Rcpp::NumericVector> size_param_x,
                bool cmHet,
                bool dmHet,
                double tolerance,
                unsigned int maxIter,
                Rcpp::Nullable<Rcpp::NumericVector> size_param_cmStart,
                Rcpp::Nullable<Rcpp::NumericVector> size_param_dmStart,
                unsigned int mG,
                unsigned int mE,
                double maxDiag,
                bool stopIfDecreasing);

void updateFAHomVar(arma::mat s,
                    arma::mat& wNew,
                    arma::mat& pNew,
                    unsigned int m,
                    double maxDiag);

void updateFA(arma::mat y,
              arma::mat wStart,
              arma::mat pStart,
              arma::mat& wNew,
              arma::mat& pNew,
              unsigned int m0,
              bool hetVar,
              double maxDiag,
              double tolerance,
              unsigned int maxIter,
              bool printProgress);

void updatePrec(unsigned int m,
                unsigned int nc,
                arma::mat omega,
                const arma::mat w,
                const arma::mat p,
                arma::mat& wNew,
                arma::mat& pNew,
                arma::mat& cNew,
                bool het,
                double maxDiag);

arma::mat vecInvDiag(arma::vec x,
                     arma::vec y);

arma::vec tracePInvDiag(arma::vec x,
                        arma::vec y);

#endif
