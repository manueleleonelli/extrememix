// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
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


// [[Rcpp::export]]
double c_rmgamma(NumericVector mu, NumericVector eta, NumericVector w){
  NumericVector cw = cumsum(w);
  int n = mu.length();
  int index = 0;
  double u = R::runif(0,1);
  if(u <= cw[0]){ index = 0;}
  for(int i = 0; i < (n-1); i++ ){
    if(u > cw[i] && u <= cw[i+1]){
      index = i+1;
    }
  }
  return(c_rgamma(mu[index],eta[index]));
}


// [[Rcpp::export]]
NumericVector c_dggpd(NumericVector x, double xi, double sigma, double u, double mu, double eta) {
  int n = x.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++){
    if(x[i] <= u){
      res[i] = c_dgamma(x[i], mu, eta);
    }
    else{ res[i] = (1-c_pgamma(u,mu,eta))*c_dgpd(x[i],xi,sigma,u);}
  }
  return res;
}

// [[Rcpp::export]]
NumericVector c_pggpd(NumericVector x, double xi, double sigma, double u, double mu, double eta) {
  int n = x.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++){
    if(x[i] <= u){
      res[i] = c_pgamma(x[i], mu, eta);
    }
    else{ res[i] = c_pgamma(u,mu,eta) + (1-c_pgamma(u,mu,eta)) * c_pgpd(x[i],xi,sigma,u);}
  }
  return res;
}


// [[Rcpp::export]]
NumericVector c_rggpd(double N, double xi, double sigma, double u, double mu, double eta){
  NumericVector res(N);
  for (int i = 0; i < N; i++){
    res[i] = c_rgamma(mu,eta);
    if(res[i]>u){res[i] = c_rgpd(xi,sigma,u);}
  }
  return res;
}

// [[Rcpp::export]]
NumericVector c_qggpd(NumericVector p, double xi, double sigma, double u, double mu, double eta) {
  int n = p.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++){
    if(p[i] <= c_pgamma(u,mu,eta)){
      res[i] = c_qgamma(p[i], mu, eta);
    }
    else{ res[i] = u + (sigma/xi)*(pow(1-(p[i]-c_pgamma(u,mu,eta))/(1-c_pgamma(u,mu,eta)),-xi)-1);}
  }
  return res;
}


// [[Rcpp::export]]
NumericVector c_dmgpd(NumericVector x, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w){
  int n = x.length();
  NumericVector out(n);
  for(int i =0; i <n; i++){
    if(x[i] <= u){out[i] = c_dmgamma(x[i],mu,eta,w);}
    else{out[i] = (1-c_pmgamma(u,mu,eta,w))*c_dgpd(x[i],xi,sigma,u);}
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector c_pmgpd(NumericVector q, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w){
  int n = q.length();
  NumericVector out(n);
  for(int i =0; i <n; i++){
    if(q[i] <= u){out[i] = c_pmgamma(q[i],mu,eta,w);}
    else{out[i] = c_pmgamma(u,mu,eta,w) + (1-c_pmgamma(u,mu,eta,w))*c_pgpd(q[i],xi,sigma,u);}
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector c_qmgpd(NumericVector p, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w){
  int n = p.length();
  NumericVector out(n);
  double thresh = c_pmgamma(u,mu,eta,w);
  for(int i =0; i <n; i++){
    if(p[i] <= thresh){out[i] = c_qmgamma(p[i],mu,eta,w);}
    else{out[i] = u + (sigma/xi)*(pow(1-(p[i]-c_pmgamma(u,mu,eta,w))/(1-c_pmgamma(u,mu,eta,w)),-xi)-1);}
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector c_rmgpd(int N, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w){
  NumericVector out(N);
  for (int i = 0; i < N; i++){
    out[i] = c_rmgamma(mu,eta,w);
    if (out[i] > u) {out[i] = c_rgpd(xi,sigma,u);}
  }
  return(out);
}

double ggpd_logprior(double xi, double sigma, double u, double mu, double eta, double mu_u, double sd_u, double a, double b, double c, double d){
  return(-log(sigma)-log(1+xi)-0.5*log(1+2*xi)-0.5*pow(u-mu_u,2)/pow(sd_u,2)+(b-1)*log(mu)-(b/a)*mu+(d-1)*log(eta)-(d/c)*eta);
}

double ggpd_logpost(NumericVector x, NumericVector param, NumericVector prior){
  double post = sum(log(c_dggpd(x,param(0),param(1),param(2),param(3),param(4))));
  return(post+ggpd_logprior(param(0),param(1),param(2),param(3),param(4),prior(0),prior(1),prior(2),prior(3),prior(4),prior(5)));
}


// [[Rcpp::export]]
List c_fggpd(NumericVector x, int it, NumericVector start, NumericVector var, NumericVector prior){
  double value = 0;
  double M = max(x);
  double m = min(x);
  double as;

  double vl;
  int l =0;

  int countxi=0;
  int countsigma=0;
  int countu=0;
  int countmu=0;
  int counteta=0;

  NumericMatrix out(it+1,5);
  NumericMatrix V(it+1,5);
  NumericVector temp(it);
  out(0,_) = start;
  V(0,_) = var;
  NumericVector cur = start;

  Progress p(it);

  for(int i = 1; i <= it; i++){

    if (Progress::check_abort() )
      stop("User Interrupted") ;
    p.increment();

    out(i,_) = out(i-1,_);
    cur = out(i,_);
    V(i,_) = V(i-1,_);

    cur(0) = out(i,0)+sqrt(V(i,0))*R::rnorm(0,1);
    while(cur(0) < -cur(1)*(M-cur(2))){cur(0) = out(i,0)+sqrt(V(i,0))*R::rnorm(0,1);}
    value = exp( ggpd_logpost(x,cur,prior) -  ggpd_logpost(x,out(i,_),prior) + R::pnorm((out(i,0)+out(i,1)/(M-out(i,2)))/sqrt(V(i,0)),0,1,1,0) - R::pnorm((cur(0)+out(i,1)/(M-out(i,2)))/sqrt(V(i,0)),0,1,1,0));
    if (value > R::runif(0,1)) {out(i,0) = cur(0); countxi += 1;} else {cur(0) = out(i-1,0);}

    if(cur(0) <0 ){
      cur(1) = out(i,1)+sqrt(V(i,1))*R::rnorm(0,1);
      while(cur(1) < -cur(0)*(M-cur(2))){cur(1) = out(i,1)+sqrt(V(i,1))*R::rnorm(0,1);}
      value = exp( ggpd_logpost(x,cur,prior) -  ggpd_logpost(x,out(i,_),prior) + R::pnorm((out(i,1)+cur(0)*(M-cur(2)))/sqrt(V(i,1)),0,1,1,0)  - R::pnorm((cur(1)+cur(0)*(M-cur(2)))/sqrt(V(i,1)),0,1,1,0) );
      if (value > R::runif(0,1)) {out(i,1) = cur(1);countsigma += 1;} else {cur(1) = out(i-1,1);}
    } else{
      cur(1) = c_rgamma(out(i,1),pow(out(i,1),2)/V(i,1));
      value = exp( ggpd_logpost(x,cur,prior) -  ggpd_logpost(x,out(i,_),prior) + c_dgamma(out(i,1),cur(1),pow(cur(1),2)/V(i,1)) - c_dgamma(cur(1),out(i,1),pow(out(i,1),2)/V(i,1)) );
      if (value > R::runif(0,1)) {out(i,1) = cur(1); countsigma += 1;} else {cur(1) = out(i-1,1);}
    }

    cur(2) = out(i,2)+sqrt(V(i,2))* R::rnorm(0,1);
    if (cur(0)<0) {as = M+ cur(1)/cur(0); } else {as = m;}
    while(cur(2) < as){cur(2) = out(i,2)+sqrt(V(i,2))*R::rnorm(0,1);}
    value = exp( ggpd_logpost(x,cur,prior) -  ggpd_logpost(x,out(i,_),prior) + R::pnorm((out(i,2)-as)/sqrt(V(i,2)),0,1,1,0) - R::pnorm((cur(2)-as)/sqrt(V(i,2)),0,1,1,0));
    if(value > R::runif(0,1)) {out(i,2) = cur(2); countu += 1;} else {cur(2) = out(i-1,2);}

    cur(3) = c_rgamma(out(i,3),pow(out(i,3),2)/V(i,3));
    value = exp(ggpd_logpost(x,cur,prior)- ggpd_logpost(x,out(i,_),prior) + c_dgamma(out(i,3),cur(3),pow(cur(3),2)/V(i,3)) - c_dgamma(cur(3),out(i,3),pow(out(i,3),2)/V(i,3)));
    if(value > R::runif(0,1)) {out(i,3) = cur(3); countmu += 1; } else {cur(3) = out(i-1,3); }

    cur(4) = c_rgamma(out(i,4),pow(out(i,4),2)/V(i,4));
    value = exp( ggpd_logpost(x,cur,prior)- ggpd_logpost(x,out(i,_),prior) +  c_dgamma(out(i,4),cur(4),pow(cur(4),2)/V(i,4))  - c_dgamma(cur(4),out(i,4),pow(out(i,4),2)/V(i,4)));
    if(value > R::runif(0,1)) {out(i,4) = cur(4); counteta += 1; } else {cur(4) = out(i-1,4); }

    if (i%50 == 0){
      if ( 0.01 <= 1/sqrt(l+1)) {vl = 0.01;} else{vl = 1/sqrt(l+1);}
      if(countxi <= 22){V(i,0) = exp(log(V(i,0)) - vl); } else {V(i,0) = exp(log(V(i,0)) + vl );}
      if(countsigma <= 22){V(i,1) = exp(log(V(i,1)) - vl); } else {V(i,1) = exp(log(V(i,1)) + vl );}
      if(countu <= 22){V(i,2) = exp(log(V(i,2)) - vl); } else {V(i,2) = exp(log(V(i,2)) + vl );}
      if(countmu <= 22){V(i,3) = exp(log(V(i,3)) - vl); } else {V(i,3) = exp(log(V(i,3)) + vl );}
      if(counteta <= 22){V(i,4) = exp(log(V(i,4)) - vl); } else {V(i,4) = exp(log(V(i,4)) + vl );}
      l += 1;
      countxi = 0;
      countsigma =0;
      countu = 0;
      countmu = 0;
      counteta = 0;
    }
  }


  List donno;
  donno["chain"] =out ;
  donno["var"] = V;
  return donno;
}


double c_mgpd_loglik(NumericVector x, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w){
  double logLik = sum(log(c_dmgpd(x,xi,sigma,u,mu,eta,w)));
  return(logLik);
}


double c_mgpd_logprior(double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector prior_u, NumericVector prior_mu, NumericVector prior_eta){
  double out = -log(sigma) - log(1+xi) -0.5 * log(1+2*xi) + -0.5*pow(((u-prior_u[0])/prior_u[1]),2);
  for (int i = 1; i < prior_mu.length(); i = i+2){
    out += log(mu[(i-1)/2])*(prior_mu[i]-1) - (prior_mu[i]/prior_mu[i-1])*mu[(i-1)/2] + log(eta[(i-1)/2])*(prior_eta[i]-1) - (prior_eta[i]/prior_eta[i-1])*eta[(i-1)/2];
  }
  return(out);
}


double c_mgpd_logpost(NumericVector x, double xi, double sigma, double u, NumericVector mu, NumericVector eta, NumericVector w, NumericVector prior_u, NumericVector prior_mu,NumericVector prior_eta){
  double out = c_mgpd_loglik(x,xi,sigma,u,mu,eta,w) + c_mgpd_logprior(xi,sigma,u, mu, eta, prior_u, prior_mu,prior_eta);
  return(out);
}

NumericVector r_dir(NumericVector w, double V){
  int k = w.length();
  NumericVector out(k);
  for(int i =0; i <k; i++){
    out[i] = c_rgamma(w[i]*V,1);
    if(out[i] < 0.00001) {out[i] = 0.00001;}
  }
  out = out/sum(out);
  return(out);
}

// [[Rcpp::export]]
List c_fmgpd(NumericVector x, int it, int k, NumericVector start_gpd, NumericVector start_mu, NumericVector start_eta, NumericVector start_w, NumericVector var, NumericVector prior_u, NumericVector prior_mu, NumericVector prior_eta){
  double value = 0;
  double M = max(x);
  double m = min(x);
  double as;

  double vl;
  int l = 0;

  int countxi=0;
  int countsigma=0;
  int countu=0;
  NumericVector countmu (k);
  int countw = 0;
  NumericVector xic(it+1); NumericVector sigmac(it+1); NumericVector uc(it+1); NumericMatrix muc(it+1,k); NumericMatrix etac(it+1,k); NumericMatrix wc(it+1,k);

  NumericMatrix out(it+1,3 +3*k); NumericMatrix V(it+1, 4 + k);

  xic[0] = start_gpd[0]; sigmac[0] = start_gpd[1]; uc[0]= start_gpd[2]; muc(0,_) = start_mu; etac(0,_) = start_eta; wc(0,_) = start_w;

  V(0,_) = var;

  double xip; double sigmap; double up; NumericVector wp; NumericVector mup; NumericVector etap;

  Progress p(it);

  for(int i = 1; i <= it; i++){

    if (Progress::check_abort() )
      stop("User Interrupted") ;
    p.increment();

    xic[i] = xic[i-1]; sigmac[i] = sigmac[i-1]; uc[i]= uc[i-1]; muc(i,_) = muc(i-1,_); etac(i,_) = etac(i-1,_); wc(i,_) = wc(i-1,_);
    V(i,_) = V(i-1,_);

    xip = xic[i]+sqrt(V(i,0))*R::rnorm(0,1);
    while(xip < -sigmac[i]*(M-uc[i])){ xip = xic[i]+sqrt(V(i,0))*R::rnorm(0,1);}
    value = exp(c_mgpd_logpost(x, xip, sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta)  + R::pnorm((xic[i]+sigmac[i]/(M-uc[i]))/sqrt(V(i,0)),0,1,1,0) - R::pnorm((xip+sigmac[i]/(M-uc[i]))/sqrt(V(i,0)),0,1,1,0));
    if (value > R::runif(0,1)) {xic[i] = xip; countxi += 1;}

    if(xic[i] <0 ){
      sigmap = sigmac[i]+sqrt(V(i,1))*R::rnorm(0,1);
      while(sigmap < -xic[i]*(M-uc[i])){sigmap = sigmac[i]+sqrt(V(i,1))*R::rnorm(0,1);}
      value = exp(c_mgpd_logpost(x, xic[i], sigmap, uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + R::pnorm((sigmac[i]+xic[i]*(M-uc[i]))/sqrt(V(i,1)),0,1,1,0)  - R::pnorm((sigmap+xic[i]*(M-uc[i]))/sqrt(V(i,1)),0,1,1,0) );
      if (value > R::runif(0,1)) {sigmac[i] = sigmap;countsigma += 1;}
    } else{
      sigmap = c_rgamma(sigmac[i],pow(sigmac[i],2)/V(i,1));
      value = exp(c_mgpd_logpost(x, xic[i], sigmap, uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + c_dgamma(sigmac[i],sigmap,pow(sigmap,2)/V(i,1)) - c_dgamma(sigmap,sigmac[i],pow(sigmac[i],2)/V(i,1)) );
      if (value > R::runif(0,1)) {sigmac[i] = sigmap;countsigma += 1;}
    }

    up = uc[i]+sqrt(V(i,2))*R::rnorm(0,1);
    if (xic[i]<0) {as = M+ sigmac[i]/xic[i]; } else {as = m;}
    while(up < as){up = uc[i]+sqrt(V(i,2))*R::rnorm(0,1);}
    value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], up, muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + R::pnorm((uc[i]-as)/sqrt(V(i,2)),0,1,1,0) - R::pnorm((up-as)/sqrt(V(i,2)),0,1,1,0));
    if (value > R::runif(0,1)) {uc[i] = up;countu += 1;}

    mup = muc(i,_); etap = etac(i,_);
    mup[0] = c_rgamma(muc(i,0),pow(muc(i,0),2)/V(i,3));
    while(mup[0] >= muc(i,1)) {mup[0] = c_rgamma(muc(i,0),pow(muc(i,0),2)/V(i,3));}
    if(mup[0] < 0.00001) {mup[0] = 0.00001;}
    etap[0] = c_rgamma(etac(i,0),pow(etac(i,0),2)/V(i,3));
    if(etap[0] < 0.00001) {etap[0] = 0.00001;}
    value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], mup, etap, wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + c_dgamma(muc(i,0),mup[0],pow(mup[0],2)/V(i,3)) -c_dgamma(mup[0],muc(i,0),pow(muc(i,0),2)/V(i,3)) + c_dgamma(etac(i,0),etap[0],pow(etap[0],2)/V(i,3)) -c_dgamma(etap[0],etac(i,0),pow(etac(i,0),2)/V(i,3)) );
    if (value > R::runif(0,1)) {muc(i,_) = mup; etac(i,_) = etap; countmu[0] += 1;} else{mup = muc(i,_); etap = etac(i,_);}


    mup[1] = c_rgamma(muc(i,1),pow(muc(i,1),2)/V(i,4));
    while(mup[1] <= muc(i,0)) {mup[1] = c_rgamma(muc(i,1),pow(muc(i,1),2)/V(i,4));}
    if(mup[1] < 0.00001) {mup[1] = 0.00001;}
    if(k >= 3) { while(mup[1]>= muc(i,2)) {mup[1] = c_rgamma(muc(i,1),pow(muc(i,1),2)/V(i,4));}}
    etap[1] = c_rgamma(etac(i,1),pow(etac(i,1),2)/V(i,4));
    if(etap[1] < 0.00001) {etap[1] = 0.00001;}
    value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], mup, etap, wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + c_dgamma(muc(i,1),mup[1],pow(mup[1],2)/V(i,4)) -c_dgamma(mup[1],muc(i,1),pow(muc(i,1),2)/V(i,4)) + c_dgamma(etac(i,1),etap[1],pow(etap[1],2)/V(i,4)) -c_dgamma(etap[1],etac(i,1),pow(etac(i,1),2)/V(i,4)) );
    if (value > R::runif(0,1)) {muc(i,_) = mup; etac(i,_) = etap; countmu[1] += 1;} else{mup = muc(i,_); etap = etac(i,_);}


    if(k >= 3){
      mup[2] = c_rgamma(muc(i,2),pow(muc(i,2),2)/V(i,5));
      while(mup[2] <= muc(i,1)) {mup[2] = c_rgamma(muc(i,2),pow(muc(i,2),2)/V(i,5));}
      if(k == 4) { while(mup[2]>= muc(i,3)) {mup[2] = c_rgamma(muc(i,2),pow(muc(i,2),2)/V(i,5));}}
      etap[2] = c_rgamma(etac(i,2),pow(etac(i,2),2)/V(i,5));
      if(etap[2] < 0.00001) {etap[2] = 0.00001;}
      value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], mup, etap, wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + c_dgamma(muc(i,2),mup[2],pow(mup[2],2)/V(i,5)) -c_dgamma(mup[2],muc(i,2),pow(muc(i,2),2)/V(i,5)) + c_dgamma(etac(i,2),etap[2],pow(etap[2],2)/V(i,5)) -c_dgamma(etap[2],etac(i,2),pow(etac(i,2),2)/V(i,5)) );
      if (value > R::runif(0,1)) {muc(i,_) = mup; etac(i,_) = etap; countmu[2] += 1;} else{mup = muc(i,_); etap = etac(i,_);}
    }

    if(k == 4){
      mup[3] = c_rgamma(muc(i,3),pow(muc(i,3),3)/V(i,6));
      while(mup[3] <= muc(i,2)) {mup[3] = c_rgamma(muc(i,3),pow(muc(i,3),2)/V(i,6));}
      etap[3] = c_rgamma(etac(i,3),pow(etac(i,3),2)/V(i,6));
      if(etap[3] < 0.00001) {etap[3] = 0.00001;}
      value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], mup, etap, wc(i,_), prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta) + c_dgamma(muc(i,3),mup[3],pow(mup[3],2)/V(i,6)) -c_dgamma(mup[3],muc(i,3),pow(muc(i,3),2)/V(i,6)) + c_dgamma(etac(i,3),etap[3],pow(etap[3],2)/V(i,6)) -c_dgamma(etap[3],etac(i,3),pow(etac(i,3),2)/V(i,6)) );
      if (value > R::runif(0,1)) {muc(i,_) = mup; etac(i,_) = etap; countmu[3] += 1;} else{mup = muc(i,_); etap = etac(i,_);}
    }



    wp = r_dir(wc(i,_),V(i,4 + k-1));
    value = exp(c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wp, prior_u, prior_mu, prior_eta) - c_mgpd_logpost(x, xic[i], sigmac[i], uc[i], muc(i,_), etac(i,_), wc(i,_), prior_u, prior_mu, prior_eta)  + sum(log(wc(i,_))) - sum(log(wp)));
    if (value > R::runif(0,1)) {wc(i,_) = wp;countw += 1;}

    if (i%50 == 0){
      if ( 0.01 <= 1/sqrt(l+1)) {vl = 0.01;} else{vl = 1/sqrt(l+1);}
      if(countxi <= 22){V(i,0) = exp(log(V(i,0)) - vl); } else {V(i,0) = exp(log(V(i,0)) + vl );}
      if(countsigma <= 22){V(i,1) = exp(log(V(i,1)) - vl); } else {V(i,1) = exp(log(V(i,1)) + vl );}
      if(countu <= 22){V(i,2) = exp(log(V(i,2)) - vl); } else {V(i,2) = exp(log(V(i,2)) + vl );}
      if(countw <= 22){V(i,3+k) = exp(log(V(i,3+k)) - vl); } else {V(i,3+k) = exp(log(V(i,3+k)) + vl );}
      if(countmu[0] <= 22){V(i,3) = exp(log(V(i,3)) - vl); } else {V(i,3) = exp(log(V(i,3)) + vl );}
      if(countmu[1] <= 22){V(i,4) = exp(log(V(i,4)) - vl); } else {V(i,4) = exp(log(V(i,4)) + vl );}
      if(k >= 3){if(countmu[2] <= 22){V(i,5) = exp(log(V(i,5)) - vl); } else {V(i,5) = exp(log(V(i,5)) + vl );}}
      if(k == 4){if(countmu[3] <= 22){V(i,6) = exp(log(V(i,6)) - vl); } else {V(i,6) = exp(log(V(i,6)) + vl );}}

      l += 1;
      countxi = 0;
      countsigma =0;
      countu = 0;
      countw = 0;
      countmu[0] = 0;
      countmu[1] = 0;
      if(k >= 3){countmu[2] = 0;}
      if(k == 4){countmu[3] = 0;}
    }

  }

  out = cbind(xic,sigmac,uc,muc,etac,wc);

  List donno;
  donno["chain"] =out ;
  donno["var"] = V;
  return donno;
}



// [[Rcpp::export]]
NumericMatrix c_quant_ggpd(NumericMatrix chain, NumericVector x){
  NumericMatrix out(chain.nrow(),x.length());
  for(int i=0; i < chain.nrow(); i++){
    out(i,_) = c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
  }
  return(out);
}

// [[Rcpp::export]]
NumericMatrix c_quant_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x){

  NumericMatrix out(mu.nrow(),x.length());
  for(int i=0; i < mu.nrow(); i++){
    out(i,_) = c_qmgpd(x,gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_), w(i,_));
  }
  return(out);
}



// [[Rcpp::export]]
double DIC_ggpd(NumericMatrix chain, NumericVector data){
  double xi = mean(chain(_,0));
  double sigma = mean(chain(_,1));
  double u = mean(chain(_,2));
  double mu = mean(chain(_,3));
  double eta = mean(chain(_,4));
  double d = -2*sum(log(c_dggpd(data,xi,sigma,u,mu,eta)));
  double dbar = 0;
  for(int i = 0; i < chain.nrow(); i++){
    dbar += -2*sum(log(c_dggpd(data,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4))));
  }
  return(2*dbar/chain.nrow() - d );
}

// [[Rcpp::export]]
double DIC_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector data){
  double xi = mean(gpd(_,0));
  double sigma = mean(gpd(_,1));
  double u = mean(gpd(_,2));
  int k = mu.ncol();
  NumericVector mup(k); NumericVector etap(k); NumericVector wp(k);
  for(int i = 0; i < k; i++){
    mup[i] = mean(mu(_,i));
    etap[i] = mean(eta(_,i));
    wp[i] = mean(w(_,i));
  }
  double d = -2*sum(log(c_dmgpd(data,xi,sigma,u,mup,etap,wp)));
  double dbar = 0;
  for(int i = 0; i < gpd.nrow(); i++){
    dbar += -2*sum(log(c_dmgpd(data,gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_),w(i,_))));
  }
  return(2*dbar/gpd.nrow() - d );
}

// [[Rcpp::export]]
double WAIC_ggpd(NumericMatrix chain, NumericVector data){
  int r = chain.nrow();
  int n = data.length();
  NumericMatrix l(r,n);
  for(int i = 0; i < r; i++){
    l(i,_) = c_dggpd(data, chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
  }
  NumericVector temp(n); NumericVector temp1(n);
  for(int i = 0; i < n; i++){temp[i] = log(sum(l(_,i))/r); temp1[i] = sum(pow(log(l(_,i))-sum(log(l(_,i)))/r,2))/(r-1);}
  return(-2*sum(temp) +2*sum(temp1) );
}

// [[Rcpp::export]]
double WAIC_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector data){
  int r = gpd.nrow();
  int n = data.length();
  NumericMatrix l(r,n);
  for(int i = 0; i < r; i++){
    l(i,_) = c_dmgpd(data, gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_),w(i,_));
  }
  NumericVector temp(n); NumericVector temp1(n);
  for(int i = 0; i < n; i++){temp[i] = log(sum(l(_,i))/r); temp1[i] = sum(pow(log(l(_,i))-sum(log(l(_,i)))/r,2))/(r-1);}
  return(-2*sum(temp) +2*sum(temp1) );
}


// [[Rcpp::export]]
NumericMatrix c_pred_ggpd(NumericVector x, NumericMatrix chain){
  NumericMatrix out(chain.nrow(),x.length());
  for(int i=0; i < chain.nrow(); i++){
    out(i,_) = c_dggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
  }
  return(out);
}


// [[Rcpp::export]]
NumericMatrix c_es_ggpd(NumericMatrix chain, NumericVector x){
  NumericMatrix out(chain.nrow(),x.length());
  for(int i=0; i < chain.nrow(); i++){
    out(i,_) = c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4))/(1-chain(i,0)) + (chain(i,1)-chain(i,0)*chain(i,2))/(1-chain(i,0));
    //out(i,_) = (chain(i,1)+chain(i,0)*c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4)))/(1-chain(i,0));
  }
  return(out);
}

// [[Rcpp::export]]
NumericMatrix c_tvar_ggpd(NumericMatrix chain, NumericVector x){
  NumericMatrix out(chain.nrow(),x.length());
  for(int i=0; i < chain.nrow(); i++){
    out(i,_) = c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4))/(1-chain(i,0)) + (chain(i,1)-chain(i,0)*chain(i,2))/(1-chain(i,0)) + c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
    //out(i,_) = (chain(i,1)+chain(i,0)*c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4)))/(1-chain(i,0))+ c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
  }
  return(out);
}

// [[Rcpp::export]]
NumericMatrix c_pred_mgpd(NumericVector x, NumericVector xi, NumericVector sigma, NumericVector u, NumericMatrix mu, NumericMatrix eta, NumericMatrix w){
  NumericMatrix out(xi.length(),x.length());
  for(int i=0; i < xi.length(); i++){
    out(i,_) = c_dmgpd(x,xi[i],sigma[i],u[i],mu(i,_),eta(i,_),w(i,_));
  }
  return(out);
}


// [[Rcpp::export]]
NumericMatrix c_es_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x){
  NumericMatrix out(mu.nrow(),x.length());
  for(int i=0; i < mu.nrow(); i++){
    out(i,_) = c_qmgpd(x,gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_), w(i,_))/(1-gpd(i,0)) + (gpd(i,1)-gpd(i,0)*gpd(i,2))/(1-gpd(i,0));
    //out(i,_) = (chain(i,1)+chain(i,0)*c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4)))/(1-chain(i,0));
  }
  return(out);
}

// [[Rcpp::export]]
NumericMatrix c_tvar_mgpd(NumericMatrix gpd, NumericMatrix mu, NumericMatrix eta, NumericMatrix w, NumericVector x){
  NumericMatrix out(mu.nrow(),x.length());
  for(int i=0; i < mu.nrow(); i++){
    out(i,_) = c_qmgpd(x,gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_), w(i,_))/(1-gpd(i,0)) + (gpd(i,1)-gpd(i,0)*gpd(i,2))/(1-gpd(i,0)) + c_qmgpd(x,gpd(i,0),gpd(i,1),gpd(i,2),mu(i,_),eta(i,_),w(i,_));
    //out(i,_) = (chain(i,1)+chain(i,0)*c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4)))/(1-chain(i,0))+ c_qggpd(x,chain(i,0),chain(i,1),chain(i,2),chain(i,3),chain(i,4));
  }
  return(out);
}



