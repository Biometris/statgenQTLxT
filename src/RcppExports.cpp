// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// updateFAHomVarCPP
void updateFAHomVarCPP(arma::mat s, arma::mat& wNew, arma::mat& pNew, unsigned int m, double maxDiag);
RcppExport SEXP _genStatPipeline_updateFAHomVarCPP(SEXP sSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP mSEXP, SEXP maxDiagSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type wNew(wNewSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pNew(pNewSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type maxDiag(maxDiagSEXP);
    updateFAHomVarCPP(s, wNew, pNew, m, maxDiag);
    return R_NilValue;
END_RCPP
}
// updateFACPP
void updateFACPP(arma::mat y, arma::mat wStart, arma::mat pStart, arma::mat& wNew, arma::mat& pNew, unsigned int m0, bool hetVar, double maxDiag, double tolerance, unsigned int maxIter, bool printProgress);
RcppExport SEXP _genStatPipeline_updateFACPP(SEXP ySEXP, SEXP wStartSEXP, SEXP pStartSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP m0SEXP, SEXP hetVarSEXP, SEXP maxDiagSEXP, SEXP toleranceSEXP, SEXP maxIterSEXP, SEXP printProgressSEXP) {
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
    updateFACPP(y, wStart, pStart, wNew, pNew, m0, hetVar, maxDiag, tolerance, maxIter, printProgress);
    return R_NilValue;
END_RCPP
}
// updatePrecCPP
void updatePrecCPP(unsigned int m, unsigned int nc, arma::mat omega, const arma::mat w, const arma::mat p, arma::mat& wNew, arma::mat& pNew, arma::mat& cNew, bool het, double maxDiag);
RcppExport SEXP _genStatPipeline_updatePrecCPP(SEXP mSEXP, SEXP ncSEXP, SEXP omegaSEXP, SEXP wSEXP, SEXP pSEXP, SEXP wNewSEXP, SEXP pNewSEXP, SEXP cNewSEXP, SEXP hetSEXP, SEXP maxDiagSEXP) {
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
    updatePrecCPP(m, nc, omega, w, p, wNew, pNew, cNew, het, maxDiag);
    return R_NilValue;
END_RCPP
}
// EMFACPP
List EMFACPP(arma::mat y, arma::mat k, Rcpp::Nullable<Rcpp::NumericVector> size_param_x, bool cmHet, bool dmHet, double tolerance, unsigned int maxIter, Rcpp::Nullable<Rcpp::NumericVector> size_param_cmStart, Rcpp::Nullable<Rcpp::NumericVector> size_param_dmStart, unsigned int mG, unsigned int mE, double maxDiag, bool stopIfDecreasing);
RcppExport SEXP _genStatPipeline_EMFACPP(SEXP ySEXP, SEXP kSEXP, SEXP size_param_xSEXP, SEXP cmHetSEXP, SEXP dmHetSEXP, SEXP toleranceSEXP, SEXP maxIterSEXP, SEXP size_param_cmStartSEXP, SEXP size_param_dmStartSEXP, SEXP mGSEXP, SEXP mESEXP, SEXP maxDiagSEXP, SEXP stopIfDecreasingSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(EMFACPP(y, k, size_param_x, cmHet, dmHet, tolerance, maxIter, size_param_cmStart, size_param_dmStart, mG, mE, maxDiag, stopIfDecreasing));
    return rcpp_result_gen;
END_RCPP
}
// LLQuadFormDiagCPP
double LLQuadFormDiagCPP(const arma::mat& y, const arma::cube& vInv, Rcpp::Nullable<Rcpp::NumericVector> size_param);
RcppExport SEXP _genStatPipeline_LLQuadFormDiagCPP(SEXP ySEXP, SEXP vInvSEXP, SEXP size_paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type vInv(vInvSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param(size_paramSEXP);
    rcpp_result_gen = Rcpp::wrap(LLQuadFormDiagCPP(y, vInv, size_param));
    return rcpp_result_gen;
END_RCPP
}
// estEffsCPP
List estEffsCPP(arma::mat y, arma::mat w, arma::mat x, arma::mat& vg, arma::mat& ve, arma::mat& k, bool returnSe, bool estCom, int ncores);
RcppExport SEXP _genStatPipeline_estEffsCPP(SEXP ySEXP, SEXP wSEXP, SEXP xSEXP, SEXP vgSEXP, SEXP veSEXP, SEXP kSEXP, SEXP returnSeSEXP, SEXP estComSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type vg(vgSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ve(veSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type returnSe(returnSeSEXP);
    Rcpp::traits::input_parameter< bool >::type estCom(estComSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(estEffsCPP(y, w, x, vg, ve, k, returnSe, estCom, ncores));
    return rcpp_result_gen;
END_RCPP
}
// fastGLSCPP
List fastGLSCPP(const arma::mat& X, const arma::vec& y, const arma::mat& sigma, Rcpp::Nullable<Rcpp::NumericVector> size_param, int ncores);
RcppExport SEXP _genStatPipeline_fastGLSCPP(SEXP XSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP size_paramSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param(size_paramSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(fastGLSCPP(X, y, sigma, size_param, ncores));
    return rcpp_result_gen;
END_RCPP
}
// fastGLSIBDCPP
List fastGLSIBDCPP(const arma::cube& mp, const arma::vec& y, const arma::mat& sigma, unsigned int ref, Rcpp::Nullable<Rcpp::NumericVector> size_param, int ncores);
RcppExport SEXP _genStatPipeline_fastGLSIBDCPP(SEXP mpSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP refSEXP, SEXP size_paramSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type mp(mpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param(size_paramSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(fastGLSIBDCPP(mp, y, sigma, ref, size_param, ncores));
    return rcpp_result_gen;
END_RCPP
}
// astleCPP
arma::mat astleCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _genStatPipeline_astleCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(astleCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// GRMCPP
arma::mat GRMCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _genStatPipeline_GRMCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(GRMCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// IBSCPP
arma::mat IBSCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _genStatPipeline_IBSCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(IBSCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// vanRadenCPP
arma::mat vanRadenCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _genStatPipeline_vanRadenCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(vanRadenCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// multiAllKinCPP
arma::mat multiAllKinCPP(arma::cube x, arma::vec posCor, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _genStatPipeline_multiAllKinCPP(SEXP xSEXP, SEXP posCorSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type posCor(posCorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(multiAllKinCPP(x, posCor, denom));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_genStatPipeline_updateFAHomVarCPP", (DL_FUNC) &_genStatPipeline_updateFAHomVarCPP, 5},
    {"_genStatPipeline_updateFACPP", (DL_FUNC) &_genStatPipeline_updateFACPP, 11},
    {"_genStatPipeline_updatePrecCPP", (DL_FUNC) &_genStatPipeline_updatePrecCPP, 10},
    {"_genStatPipeline_EMFACPP", (DL_FUNC) &_genStatPipeline_EMFACPP, 13},
    {"_genStatPipeline_LLQuadFormDiagCPP", (DL_FUNC) &_genStatPipeline_LLQuadFormDiagCPP, 3},
    {"_genStatPipeline_estEffsCPP", (DL_FUNC) &_genStatPipeline_estEffsCPP, 9},
    {"_genStatPipeline_fastGLSCPP", (DL_FUNC) &_genStatPipeline_fastGLSCPP, 5},
    {"_genStatPipeline_fastGLSIBDCPP", (DL_FUNC) &_genStatPipeline_fastGLSIBDCPP, 6},
    {"_genStatPipeline_astleCPP", (DL_FUNC) &_genStatPipeline_astleCPP, 2},
    {"_genStatPipeline_GRMCPP", (DL_FUNC) &_genStatPipeline_GRMCPP, 2},
    {"_genStatPipeline_IBSCPP", (DL_FUNC) &_genStatPipeline_IBSCPP, 2},
    {"_genStatPipeline_vanRadenCPP", (DL_FUNC) &_genStatPipeline_vanRadenCPP, 2},
    {"_genStatPipeline_multiAllKinCPP", (DL_FUNC) &_genStatPipeline_multiAllKinCPP, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_genStatPipeline(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
