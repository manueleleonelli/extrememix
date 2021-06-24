#' The GGPD distribution
#'
#' Density, distribution function, quantile function and random generation for the GGPD distribution.
#'
#'@return A co-variation matrix of the same size of the covariance matrix of \code{CI}.
#'
#'@examples dggpd(3,0.5,2,5,3,3)
#'
#'@param x vector of quantiles.
#'@param q vector of quantiles.
#'@param p vector of probabilities.
#'@param N number of observations.
#'@param xi shape parameter of the tail GPD (scalar).
#'@param sigma scale parameter of the tail GPD (scalar).
#'@param u threshold parameter of the tail GPD (scalar).
#'@param mu mean of the gamma bulk (scalar).
#'@param eta shape of the gamma bulk (scalar).
#'@param log logical; if TRUE, probabilities p are given as log(p).
#'@param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.
#'
#' @name ggpd
NULL


#' @rdname ggpd
#'
#' @export
dggpd <- function(x,xi,sigma,u,mu,eta, log = FALSE){
  if(xi < - 0.5) warning("xi is recommended to be bigger than -0.5")
  if(xi < 0 & any(x> u-sigma/xi)) stop("x beyond the upper limit")
  if(sigma <= 0) stop("sigma should be positive")
  if(any(x<=0))stop("x must be positive")
  if(mu<=0)stop("mu must be positive")
  if(eta<=0)stop("eta must be positive")
  if(u <= 0) stop("u must be positive")
  if(log == FALSE){c_dggpd(x,xi,sigma,u,mu,eta)}
  else{log(c_dggpd(x,xi,sigma,u,mu,eta))}
}

#' @rdname ggpd
#'
#' @export
pggpd <- function(q,xi,sigma,u,mu,eta,lower.tail = TRUE){
  if(xi < - 0.5) warning("xi is recommended to be bigger than -0.5")
  if(sigma <= 0) stop("sigma should be positive")
  if(mu<=0)stop("mu must be positive")
  if(eta<=0)stop("eta must be positive")
  if(u <= 0) stop("u must be positive")
  if(q< 0 || q>1){stop("q must be between zero and one")}
  if(lower.tail == TRUE){c_pggpd(q,xi,sigma,u,mu,eta)}
  else{c_pggpd(1-q,xi,sigma,u,mu,eta)}
}

#' @rdname ggpd
#'
#' @export
qggpd <- function(p,xi,sigma,u,mu,eta,lower.tail = TRUE){
  if(xi < - 0.5) warning("xi is recommended to be bigger than -0.5")
  if(sigma <= 0) stop("sigma should be positive")
  if(mu<=0)stop("mu must be positive")
  if(eta<=0)stop("eta must be positive")
   if(u <= 0) stop("u must be positive")
  if(p< 0 || p>1){stop("p must be between zero and one")}
  if(lower.tail == TRUE){c_qggpd(p,xi,sigma,u,mu,eta)}
  else{c_qggpd(1-p,xi,sigma,u,mu,eta)}
}

#' @rdname ggpd
#'
#' @export
rggpd <- function(N,xi,sigma,u,mu,eta){
  if(N%% 1 !=0)stop("N must be integer")
  if(xi < - 0.5) warning("xi is recommended to be bigger than -0.5")
  if(sigma <= 0) stop("sigma should be positive")
  if(mu<=0)stop("mu must be positive")
  if(eta<=0)stop("eta must be positive")
   if(u <= 0) stop("u must be positive")
  c_rggpd(N,xi,sigma,u,mu,eta);
}