// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_dgpd
double c_dgpd(double x, double xi, double sigma, double u);
RcppExport SEXP _extrememix_c_dgpd(SEXP xSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dgpd(x, xi, sigma, u));
    return rcpp_result_gen;
END_RCPP
}
// c_pgpd
double c_pgpd(double q, double xi, double sigma, double u);
RcppExport SEXP _extrememix_c_pgpd(SEXP qSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pgpd(q, xi, sigma, u));
    return rcpp_result_gen;
END_RCPP
}
// c_qgpd
double c_qgpd(double p, double xi, double sigma, double u);
RcppExport SEXP _extrememix_c_qgpd(SEXP pSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qgpd(p, xi, sigma, u));
    return rcpp_result_gen;
END_RCPP
}
// c_rgpd
double c_rgpd(double xi, double sigma, double u);
RcppExport SEXP _extrememix_c_rgpd(SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rgpd(xi, sigma, u));
    return rcpp_result_gen;
END_RCPP
}
// c_dgamma
double c_dgamma(double x, double mu, double eta);
RcppExport SEXP _extrememix_c_dgamma(SEXP xSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dgamma(x, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_pgamma
double c_pgamma(double q, double mu, double eta);
RcppExport SEXP _extrememix_c_pgamma(SEXP qSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pgamma(q, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_qgamma
double c_qgamma(double p, double mu, double eta);
RcppExport SEXP _extrememix_c_qgamma(SEXP pSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qgamma(p, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_rgamma
double c_rgamma(double mu, double eta);
RcppExport SEXP _extrememix_c_rgamma(SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rgamma(mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_dmgamma
double c_dmgamma(double x, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_dmgamma(SEXP xSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dmgamma(x, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_pmgamma
double c_pmgamma(double q, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_pmgamma(SEXP qSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pmgamma(q, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_qmgamma
double c_qmgamma(double p, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_qmgamma(SEXP pSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qmgamma(p, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_rmgamma
double c_rmgamma(NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_rmgamma(SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rmgamma(mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_dmgpd
NumericVector c_dmgpd(NumericVector x, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_dmgpd(SEXP xSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dmgpd(x, xi, sigma, u, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_pmgpd
NumericVector c_pmgpd(NumericVector q, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_pmgpd(SEXP qSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pmgpd(q, xi, sigma, u, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_extrememix_c_dgpd", (DL_FUNC) &_extrememix_c_dgpd, 4},
    {"_extrememix_c_pgpd", (DL_FUNC) &_extrememix_c_pgpd, 4},
    {"_extrememix_c_qgpd", (DL_FUNC) &_extrememix_c_qgpd, 4},
    {"_extrememix_c_rgpd", (DL_FUNC) &_extrememix_c_rgpd, 3},
    {"_extrememix_c_dgamma", (DL_FUNC) &_extrememix_c_dgamma, 3},
    {"_extrememix_c_pgamma", (DL_FUNC) &_extrememix_c_pgamma, 3},
    {"_extrememix_c_qgamma", (DL_FUNC) &_extrememix_c_qgamma, 3},
    {"_extrememix_c_rgamma", (DL_FUNC) &_extrememix_c_rgamma, 2},
    {"_extrememix_c_dmgamma", (DL_FUNC) &_extrememix_c_dmgamma, 4},
    {"_extrememix_c_pmgamma", (DL_FUNC) &_extrememix_c_pmgamma, 4},
    {"_extrememix_c_qmgamma", (DL_FUNC) &_extrememix_c_qmgamma, 4},
    {"_extrememix_c_rmgamma", (DL_FUNC) &_extrememix_c_rmgamma, 3},
    {"_extrememix_c_dmgpd", (DL_FUNC) &_extrememix_c_dmgpd, 7},
    {"_extrememix_c_pmgpd", (DL_FUNC) &_extrememix_c_pmgpd, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_extrememix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
