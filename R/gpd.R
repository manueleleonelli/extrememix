#' Generalized Pareto Distribution
#' @export
dgpd <- function(x, xi, sigma, mu, log = FALSE){
  if(log == FALSE){return(c_dgpd(x, xi, sigma, mu))} 
  else(return(log(c_dgpd(x,xi,sigma,mu))))
}

#' Generalized Pareto Distribution
#' @export
pgpd <- function(q, xi, sigma, mu, lower.tail = TRUE){
  if(lower.tail == TRUE){return(c_pgpd(q, xi, sigma, mu))} 
  else(return(1-c_pgpd(q,xi,sigma,mu)))
}


#' Generalized Pareto Distribution
#' @export
qgpd <- function(q, xi, sigma, mu, lower.tail = TRUE){
  if(lower.tail == TRUE){ q <- 1-q}
  return(c_qgpd(q, xi, sigma, mu))
}

#' Generalized Pareto Distribution
#' @export
rgpd <- function(N, xi, sigma, mu){
  c_rgpd(N,xi,sigma,mu)
}

#' Generalized Pareto Distribution
#' @export
dgamma <- function(x,mu,eta){
  c_dgamma(x,mu,eta)
}

#' Generalized Pareto Distribution
#' @export
pgamma <- function(q,mu,eta){
  c_pgamma(q,mu,eta)
}

#' Generalized Pareto Distribution
#' @export
qgamma <- function(p,mu,eta){
  c_qgamma(p,mu,eta)
}

#' Generalized Pareto Distribution
#' @export
rgamma <- function(mu,eta){
  c_rgamma(mu,eta)
}

#' Generalized Pareto Distribution
#' @export
dmgamma <- function(x,mu,eta,w){
  c_dmgamma(x,mu,eta,w)
}

#' Generalized Pareto Distribution
#' @export
pmgamma <- function(q,mu,eta,w){
  c_pmgamma(q,mu,eta,w)
}

#' Generalized Pareto Distribution
#' @export
qmgamma <- function(p,mu,eta,w){
  c_qmgamma(p,mu,eta,w)
}

