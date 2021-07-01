// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

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
// c_dggpd
NumericVector c_dggpd(NumericVector x, double xi, double sigma, double u, double mu, double eta);
RcppExport SEXP _extrememix_c_dggpd(SEXP xSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_dggpd(x, xi, sigma, u, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_pggpd
NumericVector c_pggpd(NumericVector x, double xi, double sigma, double u, double mu, double eta);
RcppExport SEXP _extrememix_c_pggpd(SEXP xSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pggpd(x, xi, sigma, u, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_rggpd
NumericVector c_rggpd(double N, double xi, double sigma, double u, double mu, double eta);
RcppExport SEXP _extrememix_c_rggpd(SEXP NSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rggpd(N, xi, sigma, u, mu, eta));
    return rcpp_result_gen;
END_RCPP
}
// c_qggpd
NumericVector c_qggpd(NumericVector p, double xi, double sigma, double u, double mu, double eta);
RcppExport SEXP _extrememix_c_qggpd(SEXP pSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qggpd(p, xi, sigma, u, mu, eta));
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
// c_qmgpd
NumericVector c_qmgpd(NumericVector p, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_qmgpd(SEXP pSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_qmgpd(p, xi, sigma, u, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_rmgpd
NumericVector c_rmgpd(int N, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w);
RcppExport SEXP _extrememix_c_rmgpd(SEXP NSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_rmgpd(N, xi, sigma, u, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_fmgpd
List c_fmgpd(NumericVector x, int it, int k, NumericVector start_gpd, NumericVector start_mu, NumericVector start_eta, NumericVector start_w, NumericVector var, NumericVector prior_u, NumericVector prior_mu, NumericVector prior_eta);
RcppExport SEXP _extrememix_c_fmgpd(SEXP xSEXP, SEXP itSEXP, SEXP kSEXP, SEXP start_gpdSEXP, SEXP start_muSEXP, SEXP start_etaSEXP, SEXP start_wSEXP, SEXP varSEXP, SEXP prior_uSEXP, SEXP prior_muSEXP, SEXP prior_etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start_gpd(start_gpdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start_mu(start_muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start_eta(start_etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start_w(start_wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var(varSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_u(prior_uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_mu(prior_muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior_eta(prior_etaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_fmgpd(x, it, k, start_gpd, start_mu, start_eta, start_w, var, prior_u, prior_mu, prior_eta));
    return rcpp_result_gen;
END_RCPP
}
// c_fggpd
List c_fggpd(NumericVector x, int it, NumericVector start, NumericVector var, NumericVector prior);
RcppExport SEXP _extrememix_c_fggpd(SEXP xSEXP, SEXP itSEXP, SEXP startSEXP, SEXP varSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type it(itSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var(varSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(c_fggpd(x, it, start, var, prior));
    return rcpp_result_gen;
END_RCPP
}
// DIC_ggpd
double DIC_ggpd(NumericMatrix chain, NumericVector data);
RcppExport SEXP _extrememix_DIC_ggpd(SEXP chainSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(DIC_ggpd(chain, data));
    return rcpp_result_gen;
END_RCPP
}
// WAIC_ggpd
double WAIC_ggpd(NumericMatrix chain, NumericVector data);
RcppExport SEXP _extrememix_WAIC_ggpd(SEXP chainSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(WAIC_ggpd(chain, data));
    return rcpp_result_gen;
END_RCPP
}
// c_pred_ggpd
NumericMatrix c_pred_ggpd(NumericVector x, NumericMatrix chain);
RcppExport SEXP _extrememix_c_pred_ggpd(SEXP xSEXP, SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pred_ggpd(x, chain));
    return rcpp_result_gen;
END_RCPP
}
// c_quant_ggpd
NumericMatrix c_quant_ggpd(NumericMatrix chain, NumericVector x);
RcppExport SEXP _extrememix_c_quant_ggpd(SEXP chainSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_quant_ggpd(chain, x));
    return rcpp_result_gen;
END_RCPP
}
// c_es_ggpd
NumericMatrix c_es_ggpd(NumericMatrix chain, NumericVector x);
RcppExport SEXP _extrememix_c_es_ggpd(SEXP chainSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_es_ggpd(chain, x));
    return rcpp_result_gen;
END_RCPP
}
// c_tvar_ggpd
NumericMatrix c_tvar_ggpd(NumericMatrix chain, NumericVector x);
RcppExport SEXP _extrememix_c_tvar_ggpd(SEXP chainSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_tvar_ggpd(chain, x));
    return rcpp_result_gen;
END_RCPP
}
// c_pred_mgpd
NumericMatrix c_pred_mgpd(NumericVector x, NumericVector xi, NumericVector sigma, NumericVector u, NumericMatrix mu, NumericMatrix eta, NumericMatrix w);
RcppExport SEXP _extrememix_c_pred_mgpd(SEXP xSEXP, SEXP xiSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(c_pred_mgpd(x, xi, sigma, u, mu, eta, w));
    return rcpp_result_gen;
END_RCPP
}
// c_quant_mgpd
NumericMatrix c_quant_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x);
RcppExport SEXP _extrememix_c_quant_mgpd(SEXP gpdSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gpd(gpdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_quant_mgpd(gpd, mu, eta, w, x));
    return rcpp_result_gen;
END_RCPP
}
// c_es_mgpd
NumericMatrix c_es_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x);
RcppExport SEXP _extrememix_c_es_mgpd(SEXP gpdSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gpd(gpdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_es_mgpd(gpd, mu, eta, w, x));
    return rcpp_result_gen;
END_RCPP
}
// c_tvar_mgpd
NumericMatrix c_tvar_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x);
RcppExport SEXP _extrememix_c_tvar_mgpd(SEXP gpdSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP wSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type gpd(gpdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(c_tvar_mgpd(gpd, mu, eta, w, x));
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
    {"_extrememix_c_dggpd", (DL_FUNC) &_extrememix_c_dggpd, 6},
    {"_extrememix_c_pggpd", (DL_FUNC) &_extrememix_c_pggpd, 6},
    {"_extrememix_c_rggpd", (DL_FUNC) &_extrememix_c_rggpd, 6},
    {"_extrememix_c_qggpd", (DL_FUNC) &_extrememix_c_qggpd, 6},
    {"_extrememix_c_dmgpd", (DL_FUNC) &_extrememix_c_dmgpd, 7},
    {"_extrememix_c_pmgpd", (DL_FUNC) &_extrememix_c_pmgpd, 7},
    {"_extrememix_c_qmgpd", (DL_FUNC) &_extrememix_c_qmgpd, 7},
    {"_extrememix_c_rmgpd", (DL_FUNC) &_extrememix_c_rmgpd, 7},
    {"_extrememix_c_fmgpd", (DL_FUNC) &_extrememix_c_fmgpd, 11},
    {"_extrememix_c_fggpd", (DL_FUNC) &_extrememix_c_fggpd, 5},
    {"_extrememix_DIC_ggpd", (DL_FUNC) &_extrememix_DIC_ggpd, 2},
    {"_extrememix_WAIC_ggpd", (DL_FUNC) &_extrememix_WAIC_ggpd, 2},
    {"_extrememix_c_pred_ggpd", (DL_FUNC) &_extrememix_c_pred_ggpd, 2},
    {"_extrememix_c_quant_ggpd", (DL_FUNC) &_extrememix_c_quant_ggpd, 2},
    {"_extrememix_c_es_ggpd", (DL_FUNC) &_extrememix_c_es_ggpd, 2},
    {"_extrememix_c_tvar_ggpd", (DL_FUNC) &_extrememix_c_tvar_ggpd, 2},
    {"_extrememix_c_pred_mgpd", (DL_FUNC) &_extrememix_c_pred_mgpd, 7},
    {"_extrememix_c_quant_mgpd", (DL_FUNC) &_extrememix_c_quant_mgpd, 5},
    {"_extrememix_c_es_mgpd", (DL_FUNC) &_extrememix_c_es_mgpd, 5},
    {"_extrememix_c_tvar_mgpd", (DL_FUNC) &_extrememix_c_tvar_mgpd, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_extrememix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
