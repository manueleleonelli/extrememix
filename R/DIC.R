#' Deviance Information Criterion
#'
#' Computation of the DIC for an extreme value mixture model
#'
#' Let \eqn{y} denote a dataset and \eqn{p(y|\theta)} the likelihood of a parametric model with parameter \eqn{\theta}. The deviance is defined as \eqn{D(\theta)= -2\log p(y|\theta)}. The deviance information criterion (DIC) is defined as \deqn{DIC = D(\hat\theta) + 2p_D,} where \eqn{\hat\theta} is the posterior estimate of \eqn{\theta} and \eqn{p_D} is referred to as the effective number of parameters and defined as \deqn{E_{\theta|y}(D(\theta)) - D(\hat\theta).} Models with a smaller DIC are favored.
#'
#' @param x the output of a model estimated with \code{extrememix}
#' @param ... additional arguments for compatibility.
#' @name DIC
#' @return The DIC of a model estimated with \code{extrememix}
#' @references Spiegelhalter, David J., et al. "Bayesian measures of model complexity and fit." Journal of the Royal Statistical Society: Series B 64.4 (2002): 583-639.
#' @examples DIC(rainfall_ggpd)
#'
#' @seealso \code{\link{WAIC}}
#'
#'@export
#'
DIC <- function (x, ...) {
  UseMethod("DIC", x)
}



#' @method DIC evmm
#'@export
#' @rdname DIC
#'
DIC.evmm <- function(x,...){
  if(sum(class(x) == "ggpd") == 1) return(DIC_ggpd(x$chain,x$data))
  if(sum(class(x) == "mgpd") == 1) {  k <- (ncol(x$chain)-3)/3
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  return(DIC_mgpd(gpd,mu,eta,w,x$data))
  }
}
