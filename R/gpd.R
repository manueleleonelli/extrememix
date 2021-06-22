
dgpd <- function(x, xi, sigma, u){
  return(log(c_dgpd(x,xi,sigma,u)))
}


pgpd <- function(q, xi, sigma, u){
return(1-c_pgpd(q,xi,sigma,u))
}



qgpd <- function(q, xi, sigma, u){
  return(c_qgpd(q, xi, sigma, u))
}


rgpd <- function(N, xi, sigma, u){
  c_rgpd(xi,sigma,u)
}


dgamma <- function(x,mu,eta){
  c_dgamma(x,mu,eta)
}


pgamma <- function(q,mu,eta){
  c_pgamma(q,mu,eta)
}


qgamma <- function(p,mu,eta){
  c_qgamma(p,mu,eta)
}


rgamma <- function(mu,eta){
  c_rgamma(mu,eta)
}


dmgamma <- function(x,mu,eta,w){
  c_dmgamma(x,mu,eta,w)
}


pmgamma <- function(q,mu,eta,w){
  c_pmgamma(q,mu,eta,w)
}


qmgamma <- function(p,mu,eta,w){
  c_qmgamma(p,mu,eta,w)
}



rmgamma <- function(mu,eta,w){
  c_rmgamma(mu,eta,w);
}


