#' The Gamma Mixture Distribution
#'
#' Density, distribution function, quantile function and random generation for the mixture of Gamma distribution.
#'
#' The Gamma distribution has density \deqn{f_{GA}(x|\mu,\eta)= \frac{(\eta/\mu)^\eta}{\Gamma(\eta)}x^{\eta-1}\exp(-(\eta/\mu)x), \hspace{1cm} x>0,} where \eqn{\mu>0} is the mean of the distribution and \eqn{\eta>0} is its shape.
#'  The density of a mixture of Gamma distributions with \eqn{k} components is defined as  \deqn{f_{MG}(x|\mu,\eta,w)=\sum_{i=1}^k w_if_{GA}(x|\mu_i,\eta_i),} where \eqn{w_i,\mu_i,\eta_i >0}, for \eqn{i=1,\dots,k}, \eqn{w_1+\cdots+w_k=1}, \eqn{\mu=(\mu_1,\dots,\mu_k)}, \eqn{\eta = (\eta_1,\dots,\eta_k)} and \eqn{w=(w_1,\dots,w_k)}.
#'
#'@return \code{dmgamma} gives the density, \code{pmgamma} gives the distribution function, \code{qmgamma} gives the quantile function, and \code{rmgamma} generates random deviates.
#'
#'The length of the result is determined by \code{N} for \code{rmgamma} and by the length of \code{x}, \code{q} or \code{p} otherwise.
#'
#'@references Wiper, Michael, David Rios Insua, and Fabrizio Ruggeri. "Mixtures of gamma distributions with applications." Journal of Computational and Graphical Statistics 10.3 (2001): 440-454.
#'
#'@examples dmgamma(3, mu = c(2,3), eta = c(1,2), w = c(0.3,0.7))
#'
#'@param x,q vector of quantiles.
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

