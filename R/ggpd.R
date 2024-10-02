#' The GGPD distribution
#'
#' Density, distribution function, quantile function and random generation for the GGPD distribution.
#'
#'@return  The GGPD distribution is an extreme value mixture model with density
#' \deqn{f_{GGPD}(x|\xi,\sigma,u,\mu,\eta,w)=\left\{\begin{array}{ll} f_{GA}(x|\mu,\eta), & x\leq u \\ (1-F_{GA}(u|\mu,\eta))f_{GPD}(x|\xi,\sigma,u), &\mbox{otherwise},  \end{array}\right.} where \eqn{f_{GA}} is the density of the Gamma parametrized by mean \eqn{\mu} and shape \eqn{\eta}, \eqn{F_{GA}} is the distribution function of the Gamma and \eqn{f_{GPD}} is the density of the Generalized Pareto Distribution, i.e.
#'  \deqn{f_{GPD}(x|\xi,\sigma,u)=\left\{\begin{array}{ll} 1- (1+\frac{\xi}{\sigma}(x-u))^{-1/\xi}, & \mbox{if } \xi\neq 0,\\ 1- \exp\left(-\frac{x-u}{\sigma}\right), & \mbox{if } \xi = 0, \end{array}\right.}
#' where \eqn{\xi} is a shape parameter, \eqn{\sigma > 0} is a scale parameter and \eqn{u>0} is a threshold.
#'
#' @return \code{dggpd} gives the density, \code{pggpd} gives the distribution function, \code{qggpd} gives the quantile function, and \code{rggpd} generates random deviates. The length of the result is determined by \code{N} for \code{rggpd} and by the length of \code{x}, \code{q} or \code{p} otherwise.
#'
#' @references Behrens, Cibele N., Hedibert F. Lopes, and Dani Gamerman. "Bayesian analysis of extreme events with threshold estimation." Statistical Modelling 4.3 (2004): 227-244.
#'
#'@examples dggpd(3, xi = 0.5, sigma = 2, u = 5, mu = 3, eta = 3)
#'
#'
#'@param x,q vector of quantiles.
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
  if(lower.tail == TRUE){c_pggpd(q,xi,sigma,u,mu,eta)}
  else{1-c_pggpd(q,xi,sigma,u,mu,eta)}
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
