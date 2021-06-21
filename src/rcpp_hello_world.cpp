#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double c_dgpd(double x, double xi, double sigma, double u){
  if ( sigma <= 0.0) {stop("'sigma' must be positive");}
  if (x< u) {stop("'u' must be smaller than 'x'");}
  if (fabs(xi) > 0.000001){ return(exp(-log(sigma)-((1+xi)/xi)*log(1+(xi/sigma)*(x-u))));}
  else { return(exp(-log(sigma)-(x-u)/sigma));}
}


// [[Rcpp::export]]
double c_pgpd(double q, double xi, double sigma, double u){
  if (sigma <= 0.0) {stop("'sigma' must be positive");}
  if (q< u) {stop("'u' must be smaller than 'q'");}
  if (fabs(xi) > 0.000001){return(1-pow(1+xi*(q-u)/sigma,-1/xi));}
  else{return(1-exp(-(q-u)/sigma));}
}

// [[Rcpp::export]]
double c_qgpd(double p, double xi, double sigma, double u){
  if (sigma <= 0.0) {stop("'sigma' must be positive");}
  if (p< 0) {stop("'p' must be bigger than zero");}
  if (p> 1) {stop("'p' must be smaller than one");}
  if (fabs(xi) > 0.000001){return(u + (sigma/xi)*(pow(p,-xi)-1));}
  else{return(u -sigma*log(p));}
}

// [[Rcpp::export]]
double c_rgpd(double xi, double sigma, double u){
  if (sigma <= 0.0){stop("'sigma' must be positive");}
  if (fabs(xi) > 0.000001){return(u + (sigma/xi)*(pow(R::runif(0,1),-xi)-1));}
  else{return(u -sigma*log(R::runif(0,1)));}
}

// [[Rcpp::export]]
double c_dgamma(double x, double mu, double eta){
  return(R::dgamma(x, eta, mu/eta, false));}

// [[Rcpp::export]]
double c_pgamma(double q, double mu, double eta){
  return(R::pgamma(q, eta, mu/eta, true, false));}

// [[Rcpp::export]]
double c_qgamma(double p, double mu, double eta){
  return(R::qgamma(p, eta, mu/eta, true, false));}

// [[Rcpp::export]]
double c_rgamma(double mu, double eta){
  return(R::rgamma(eta, mu/eta));}

// [[Rcpp::export]]
double c_dmgamma(double x, NumericVector mu, NumericVector eta, NumericVector w){
  int n = mu.length();
  double out = 0;
  for (int i = 0; i < n; i++){
    out = out +w[i]*c_dgamma(x,mu[i],eta[i]);
  }
  return(out);
  }

// [[Rcpp::export]]
double c_pmgamma(double q, NumericVector mu, NumericVector eta, NumericVector w){
  int n = mu.length();
  double out = 0;
  for (int i = 0; i < n; i++){
    out = out +w[i]*c_pgamma(q,mu[i],eta[i]);
  }
  return(out);
}

// [[Rcpp::export]]
double c_qmgamma(double p, NumericVector mu, NumericVector eta, NumericVector w){
  double a = 0.001;
  double b = c_qgamma(0.99999, max(mu), 0.01);
  int N = 1;
  double c = 0;
  do{
    c = (a+b)/2;
    double value = c_pmgamma(c,mu,eta,w);
    if (fabs(value - p) < 0.000001) {return(c);}
    N = N + 1;
    if ((c_pmgamma(a,mu,eta,w) - p > 0) == (value  - p > 0)) { a = c; } else {b = c;}
  } while (N <= 100000);
  stop("error");
}
