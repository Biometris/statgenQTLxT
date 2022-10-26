#ifndef NEARESTPD_H_
#define NEARESTPD_H_

#include "Rcpp.h"

arma::mat nearestPDQTLxT(arma::mat x, const bool corr, const bool keepDiag, const bool do2eigen, const bool doSym, const bool doDykstra, const double eigTol, const double convTol, const double posdTol, const int maxIter);

#endif
