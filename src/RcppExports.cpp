// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// updateFAHomVar
void updateFAHomVar(arma::mat s, arma::mat& wNew, arma::mat& pNew, unsigned int m, double maxDiag);
RcppExport SEXP _statgenQTLxT_updateFAHomVar(SEXP sSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP mSEXP, SEXP maxDiagSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type wNew(wNewSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pNew(pNewSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type maxDiag(maxDiagSEXP);
    updateFAHomVar(s, wNew, pNew, m, maxDiag);
    return R_NilValue;
END_RCPP
}
// updateFA
void updateFA(arma::mat y, arma::mat wStart, arma::mat pStart, arma::mat& wNew, arma::mat& pNew, unsigned int m0, bool hetVar, double maxDiag, double tolerance, unsigned int maxIter, bool printProgress);
RcppExport SEXP _statgenQTLxT_updateFA(SEXP ySEXP, SEXP wStartSEXP, SEXP pStartSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP m0SEXP, SEXP hetVarSEXP, SEXP maxDiagSEXP, SEXP toleranceSEXP, SEXP maxIterSEXP, SEXP printProgressSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type wStart(wStartSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pStart(pStartSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type wNew(wNewSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pNew(pNewSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< bool >::type hetVar(hetVarSEXP);
    Rcpp::traits::input_parameter< double >::type maxDiag(maxDiagSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type printProgress(printProgressSEXP);
    updateFA(y, wStart, pStart, wNew, pNew, m0, hetVar, maxDiag, tolerance, maxIter, printProgress);
    return R_NilValue;
END_RCPP
}
// updatePrec
void updatePrec(unsigned int m, unsigned int nc, arma::mat omega, const arma::mat w, const arma::mat p, arma::mat& wNew, arma::mat& pNew, arma::mat& cNew, bool het, double maxDiag);
RcppExport SEXP _statgenQTLxT_updatePrec(SEXP mSEXP, SEXP ncSEXP, SEXP omegaSEXP, SEXP wSEXP, SEXP pSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP cNewSEXP, SEXP hetSEXP, SEXP maxDiagSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type m(mSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type wNew(wNewSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pNew(pNewSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cNew(cNewSEXP);
    Rcpp::traits::input_parameter< bool >::type het(hetSEXP);
    Rcpp::traits::input_parameter< double >::type maxDiag(maxDiagSEXP);
    updatePrec(m, nc, omega, w, p, wNew, pNew, cNew, het, maxDiag);
    return R_NilValue;
END_RCPP
}
// vecInvDiag
arma::mat vecInvDiag(arma::vec x, arma::vec y);
RcppExport SEXP _statgenQTLxT_vecInvDiag(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vecInvDiag(x, y));
    return rcpp_result_gen;
END_RCPP
}
// tracePInvDiag
arma::vec tracePInvDiag(arma::vec x, arma::vec y);
RcppExport SEXP _statgenQTLxT_tracePInvDiag(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(tracePInvDiag(x, y));
    return rcpp_result_gen;
END_RCPP
}
// EMFA
List EMFA(arma::mat y, arma::mat k, Rcpp::Nullable<Rcpp::NumericVector> size_param_x, bool cmHet, bool dmHet, double tolerance, unsigned int maxIter, Rcpp::Nullable<Rcpp::NumericVector> size_param_cmStart, Rcpp::Nullable<Rcpp::NumericVector> size_param_dmStart, unsigned int mG, unsigned int mE, double maxDiag, bool stopIfDecreasing, Rcpp::CharacterVector traits);
RcppExport SEXP _statgenQTLxT_EMFA(SEXP ySEXP, SEXP kSEXP, SEXP size_param_xSEXP, SEXP cmHetSEXP, SEXP dmHetSEXP, SEXP toleranceSEXP, SEXP maxIterSEXP, SEXP size_param_cmStartSEXP, SEXP size_param_dmStartSEXP, SEXP mGSEXP, SEXP mESEXP, SEXP maxDiagSEXP, SEXP stopIfDecreasingSEXP, SEXP traitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param_x(size_param_xSEXP);
    Rcpp::traits::input_parameter< bool >::type cmHet(cmHetSEXP);
    Rcpp::traits::input_parameter< bool >::type dmHet(dmHetSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param_cmStart(size_param_cmStartSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param_dmStart(size_param_dmStartSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type mG(mGSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type mE(mESEXP);
    Rcpp::traits::input_parameter< double >::type maxDiag(maxDiagSEXP);
    Rcpp::traits::input_parameter< bool >::type stopIfDecreasing(stopIfDecreasingSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type traits(traitsSEXP);
    rcpp_result_gen = Rcpp::wrap(EMFA(y, k, size_param_x, cmHet, dmHet, tolerance, maxIter, size_param_cmStart, size_param_dmStart, mG, mE, maxDiag, stopIfDecreasing, traits));
    return rcpp_result_gen;
END_RCPP
}
// LLQuadFormDiagCPP
double LLQuadFormDiagCPP(const arma::mat& y, const arma::cube& vInv, Rcpp::Nullable<Rcpp::NumericVector> size_param_x);
RcppExport SEXP _statgenQTLxT_LLQuadFormDiagCPP(SEXP ySEXP, SEXP vInvSEXP, SEXP size_param_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type vInv(vInvSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param_x(size_param_xSEXP);
    rcpp_result_gen = Rcpp::wrap(LLQuadFormDiagCPP(y, vInv, size_param_x));
    return rcpp_result_gen;
END_RCPP
}
// estEffsCPP
List estEffsCPP(arma::mat& y0, arma::mat& w0, arma::mat& x0, arma::mat& vg, arma::mat& ve, arma::mat& k, bool returnSe, bool estCom, Rcpp::Nullable<Rcpp::IntegerVector> nCores);
RcppExport SEXP _statgenQTLxT_estEffsCPP(SEXP y0SEXP, SEXP w0SEXP, SEXP x0SEXP, SEXP vgSEXP, SEXP veSEXP, SEXP kSEXP, SEXP returnSeSEXP, SEXP estComSEXP, SEXP nCoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type vg(vgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ve(veSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type returnSe(returnSeSEXP);
    Rcpp::traits::input_parameter< bool >::type estCom(estComSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type nCores(nCoresSEXP);
    rcpp_result_gen = Rcpp::wrap(estEffsCPP(y0, w0, x0, vg, ve, k, returnSe, estCom, nCores));
    return rcpp_result_gen;
END_RCPP
}
// getThr
int getThr(Rcpp::Nullable<Rcpp::IntegerVector> nCores);
RcppExport SEXP _statgenQTLxT_getThr(SEXP nCoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type nCores(nCoresSEXP);
    rcpp_result_gen = Rcpp::wrap(getThr(nCores));
    return rcpp_result_gen;
END_RCPP
}
// nearestPD2
arma::mat nearestPD2(arma::mat x, const bool corr, const bool keepDiag, const bool do2eigen, const bool doSym, const bool doDykstra, const double eigTol, const double convTol, const double posdTol, const int maxIter);
RcppExport SEXP _statgenQTLxT_nearestPD2(SEXP xSEXP, SEXP corrSEXP, SEXP keepDiagSEXP, SEXP do2eigenSEXP, SEXP doSymSEXP, SEXP doDykstraSEXP, SEXP eigTolSEXP, SEXP convTolSEXP, SEXP posdTolSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const bool >::type keepDiag(keepDiagSEXP);
    Rcpp::traits::input_parameter< const bool >::type do2eigen(do2eigenSEXP);
    Rcpp::traits::input_parameter< const bool >::type doSym(doSymSEXP);
    Rcpp::traits::input_parameter< const bool >::type doDykstra(doDykstraSEXP);
    Rcpp::traits::input_parameter< const double >::type eigTol(eigTolSEXP);
    Rcpp::traits::input_parameter< const double >::type convTol(convTolSEXP);
    Rcpp::traits::input_parameter< const double >::type posdTol(posdTolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(nearestPD2(x, corr, keepDiag, do2eigen, doSym, doDykstra, eigTol, convTol, posdTol, maxIter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_statgenQTLxT_updateFAHomVar", (DL_FUNC) &_statgenQTLxT_updateFAHomVar, 5},
    {"_statgenQTLxT_updateFA", (DL_FUNC) &_statgenQTLxT_updateFA, 11},
    {"_statgenQTLxT_updatePrec", (DL_FUNC) &_statgenQTLxT_updatePrec, 10},
    {"_statgenQTLxT_vecInvDiag", (DL_FUNC) &_statgenQTLxT_vecInvDiag, 2},
    {"_statgenQTLxT_tracePInvDiag", (DL_FUNC) &_statgenQTLxT_tracePInvDiag, 2},
    {"_statgenQTLxT_EMFA", (DL_FUNC) &_statgenQTLxT_EMFA, 14},
    {"_statgenQTLxT_LLQuadFormDiagCPP", (DL_FUNC) &_statgenQTLxT_LLQuadFormDiagCPP, 3},
    {"_statgenQTLxT_estEffsCPP", (DL_FUNC) &_statgenQTLxT_estEffsCPP, 9},
    {"_statgenQTLxT_getThr", (DL_FUNC) &_statgenQTLxT_getThr, 1},
    {"_statgenQTLxT_nearestPD2", (DL_FUNC) &_statgenQTLxT_nearestPD2, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_statgenQTLxT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
