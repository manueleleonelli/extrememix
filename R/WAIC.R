#' Widely Applicable Information Criteria
#'
#' Computation of the WAIC for an extreme value mixture model.
#'
#' Consider a dataset \eqn{y=(y_1,\dots,y_n)}, \eqn{p(y|\theta)} the likelihood of a parametric model with parameter \eqn{\theta}, and \eqn{(\theta^{(1)},\dots,\theta^{(S)})} a sample from the posterior distribution \eqn{p(\theta|y)}.
#' Define \deqn{\textnormal{llpd} = \sum_{i=1}^n \log\left(\sum_{i=1}^Sp(y_i|\theta^{(s)}\right)} and \deqn{p_\textnormal{WAIC} = \sum_{i=1}^n Var_{\theta|y}(\log p(y_i|\theta)).}
#' Then the Widely Applicable Information Criteria is defined as \deqn{WAIC = -2\textnormal{llpd} + 2p_\textnormal{WAIC}.} Models with a smaller WAIC are favored.
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name WAIC
#' @seealso \code{\link{DIC}}
#' @return The WAIC of a model estimated with \code{extrememix}
#' @examples WAIC(rainfall_ggpd)
#'
#' @references Gelman, Andrew, Jessica Hwang, and Aki Vehtari. "Understanding predictive information criteria for Bayesian models." Statistics and computing 24.6 (2014): 997-1016.
#' @references Watanabe, Sumio. "A widely applicable Bayesian information criterion." Journal of Machine Learning Research 14.Mar (2013): 867-897.
#'
#' @export
WAIC <- function (x, ...) {
  UseMethod("WAIC", x)
}



#' @method WAIC evmm
#'@export
#' @rdname WAIC
#'
WAIC.evmm <- function(x,...){
  if(sum(class(x) == "ggpd") == 1)  return(WAIC_ggpd(x$chain,x$data))
  if(sum(class(x) == "mgpd") == 1) {k <- (ncol(x$chain)-3)/3
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  return(WAIC_mgpd(gpd,mu,eta,w,x$data))}
}

