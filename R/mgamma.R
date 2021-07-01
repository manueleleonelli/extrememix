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

#'@param mu means of the gamma mixture components (vector).
#'@param eta shapes of the gamma mixture components (vector).
#'@param w weights of the gamma mixture components (vector). Must sum to one.
#'@param log logical; if TRUE, probabilities p are given as log(p).
#'@param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.
#'
#' @name mgamma
NULL

#' @rdname mgamma
#'
#' @export
dmgamma <- function(x,mu,eta,w, log = FALSE){
  if(log == FALSE)sapply(x, function(x) c_dmgamma(x,mu,eta,w))
  else{log(sapply(c_dmgamma(x,mu,eta,w)))}
}

#' @rdname mgamma
#'
#' @export
pmgamma <- function(q,mu,eta,w, lower.tail = TRUE){
  if(lower.tail == TRUE) sapply(q, function(q)c_pmgamma(q,mu,eta,w))
  else{ sapply(q, function(q) c_pmgamma(1-q,mu,eta,w))}
}

#' @rdname mgamma
#'
#' @export
qmgamma <- function(p,mu,eta,w, lower.tail = TRUE){
  if(lower.tail == TRUE) sapply(p, function(p)c_qmgamma(p,mu,eta,w))
  else {sapply(p,function(p) c_qmgamma(1-p,mu,eta,w))}
}


#' @rdname mgamma
#'
#' @export
rmgamma <- function(N,mu,eta,w){
  replicate(N,c_rmgamma(mu,eta,w));
}

