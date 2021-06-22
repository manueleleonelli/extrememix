#' The MGPD distribution
#'
#' Density, distribution function, quantile function and random generation for the MGPD distribution.
#'
#'@return A co-variation matrix of the same size of the covariance matrix of \code{CI}.
#'
#'@examples dmgpd(3,0.5,2,5,c(2,3),c(1,2),c(0.3,0.7))
#'
#'@param x vector of quantiles.
#'@param q vector of quantiles.
#'@param p vector of probabilities.
#'@param N number of observations.
#'@param xi shape parameter of the tail GPD (scalar).
#'@param sigma scale parameter of the tail GPD (scalar).
#'@param u threshold parameter of the tail GPD (scalar).
#'@param mu means of the gamma mixture components (vector).
#'@param eta shapes of the gamma mixture components (vector).
#'@param w weights of the gamma mixture components (vector). Must sum to one.
#'@param log logical; if TRUE, probabilities p are given as log(p).
#'@param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.
#'
#' @name mgpd
NULL


#' @rdname mgpd
#'
#' @export
dmgpd <- function(x,xi,sigma,u,mu,eta,w, log = FALSE){
  if(log == FALSE){c_dmgpd(x,xi,sigma,u,mu,eta,w)}
  else{log(c_dmgpd(x,xi,sigma,u,mu,eta,w))}
}

#' @rdname mgpd
#'
#' @export
pmgpd <- function(q,xi,sigma,u,mu,eta,w, lower.tail = TRUE){
  if(lower.tail == TRUE){c_pmgpd(q,xi,sigma,u,mu,eta,w)}
  else{c_pmgpd(1-q,xi,sigma,u,mu,eta,w)}
}

#' @rdname mgpd
#'
#' @export
qmgpd <- function(p,xi,sigma,u,mu,eta,w,lower.tail = TRUE){
  if(lower.tail == TRUE){c_qmgpd(p,xi,sigma,u,mu,eta,w)}
  else{c_qmgpd(1-p,xi,sigma,u,mu,eta,w)}
}

#' @rdname mgpd
#'
#' @export
rmgpd <- function(N,xi,sigma,u,mu,eta,w){
  c_rmgpd(N,xi,sigma,u,mu,eta,w);
}