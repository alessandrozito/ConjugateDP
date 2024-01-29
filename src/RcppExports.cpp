// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_pdf_Sg
double log_pdf_Sg(double x, double a, double b, double m);
RcppExport SEXP _ConjugateDP_log_pdf_Sg(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(log_pdf_Sg(x, a, b, m));
    return rcpp_result_gen;
END_RCPP
}
// rSg_ratio_uniforms
arma::vec rSg_ratio_uniforms(int nsamples, double a, double b, double m, double Mu, double Mv);
RcppExport SEXP _ConjugateDP_rSg_ratio_uniforms(SEXP nsamplesSEXP, SEXP aSEXP, SEXP bSEXP, SEXP mSEXP, SEXP MuSEXP, SEXP MvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< double >::type Mv(MvSEXP);
    rcpp_result_gen = Rcpp::wrap(rSg_ratio_uniforms(nsamples, a, b, m, Mu, Mv));
    return rcpp_result_gen;
END_RCPP
}
// rbeta_cpp
double rbeta_cpp(double a0, double b0);
RcppExport SEXP _ConjugateDP_rbeta_cpp(SEXP a0SEXP, SEXP b0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    rcpp_result_gen = Rcpp::wrap(rbeta_cpp(a0, b0));
    return rcpp_result_gen;
END_RCPP
}
// rSg_beta_prime
arma::vec rSg_beta_prime(int nsamples, double a, double b, double m);
RcppExport SEXP _ConjugateDP_rSg_beta_prime(SEXP nsamplesSEXP, SEXP aSEXP, SEXP bSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rSg_beta_prime(nsamples, a, b, m));
    return rcpp_result_gen;
END_RCPP
}
// log_pdf_Sg_posterior
double log_pdf_Sg_posterior(double x, double a, double b, double m, double k, double n);
RcppExport SEXP _ConjugateDP_log_pdf_Sg_posterior(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP mSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(log_pdf_Sg_posterior(x, a, b, m, k, n));
    return rcpp_result_gen;
END_RCPP
}
// rSg_posterior_ratio_uniforms
arma::vec rSg_posterior_ratio_uniforms(int nsamples, double a, double b, double m, double k, double n, double Mu, double Mv);
RcppExport SEXP _ConjugateDP_rSg_posterior_ratio_uniforms(SEXP nsamplesSEXP, SEXP aSEXP, SEXP bSEXP, SEXP mSEXP, SEXP kSEXP, SEXP nSEXP, SEXP MuSEXP, SEXP MvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< double >::type Mv(MvSEXP);
    rcpp_result_gen = Rcpp::wrap(rSg_posterior_ratio_uniforms(nsamples, a, b, m, k, n, Mu, Mv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ConjugateDP_log_pdf_Sg", (DL_FUNC) &_ConjugateDP_log_pdf_Sg, 4},
    {"_ConjugateDP_rSg_ratio_uniforms", (DL_FUNC) &_ConjugateDP_rSg_ratio_uniforms, 6},
    {"_ConjugateDP_rbeta_cpp", (DL_FUNC) &_ConjugateDP_rbeta_cpp, 2},
    {"_ConjugateDP_rSg_beta_prime", (DL_FUNC) &_ConjugateDP_rSg_beta_prime, 4},
    {"_ConjugateDP_log_pdf_Sg_posterior", (DL_FUNC) &_ConjugateDP_log_pdf_Sg_posterior, 6},
    {"_ConjugateDP_rSg_posterior_ratio_uniforms", (DL_FUNC) &_ConjugateDP_rSg_posterior_ratio_uniforms, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_ConjugateDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
