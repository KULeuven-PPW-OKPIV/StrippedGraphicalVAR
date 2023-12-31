// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// beta_ridge_C
NumericMatrix beta_ridge_C(NumericMatrix X, NumericMatrix Y, double lambda_beta);
RcppExport SEXP _StrippedGraphicalVAR_beta_ridge_C(SEXP XSEXP, SEXP YSEXP, SEXP lambda_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_beta(lambda_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_ridge_C(X, Y, lambda_beta));
    return rcpp_result_gen;
END_RCPP
}
// Beta_C
NumericMatrix Beta_C(NumericMatrix kappa, NumericMatrix beta, NumericMatrix X, NumericMatrix Y, double lambda_beta, NumericMatrix lambda_beta_mat, double convergence, int maxit);
RcppExport SEXP _StrippedGraphicalVAR_Beta_C(SEXP kappaSEXP, SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP lambda_betaSEXP, SEXP lambda_beta_matSEXP, SEXP convergenceSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_beta(lambda_betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda_beta_mat(lambda_beta_matSEXP);
    Rcpp::traits::input_parameter< double >::type convergence(convergenceSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(Beta_C(kappa, beta, X, Y, lambda_beta, lambda_beta_mat, convergence, maxit));
    return rcpp_result_gen;
END_RCPP
}
// inner_loop
double inner_loop(const NumericMatrix& S, NumericMatrix& W, NumericMatrix& X, NumericVector& WXj, const NumericVector& Wd, int n, NumericMatrix& rho, double thrLasso);
RcppExport SEXP _StrippedGraphicalVAR_inner_loop(SEXP SSEXP, SEXP WSEXP, SEXP XSEXP, SEXP WXjSEXP, SEXP WdSEXP, SEXP nSEXP, SEXP rhoSEXP, SEXP thrLassoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type W(WSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type WXj(WXjSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Wd(WdSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type thrLasso(thrLassoSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_loop(S, W, X, WXj, Wd, n, rho, thrLasso));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StrippedGraphicalVAR_beta_ridge_C", (DL_FUNC) &_StrippedGraphicalVAR_beta_ridge_C, 3},
    {"_StrippedGraphicalVAR_Beta_C", (DL_FUNC) &_StrippedGraphicalVAR_Beta_C, 8},
    {"_StrippedGraphicalVAR_inner_loop", (DL_FUNC) &_StrippedGraphicalVAR_inner_loop, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_StrippedGraphicalVAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
