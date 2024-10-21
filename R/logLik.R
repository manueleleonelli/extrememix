#' Log-likelihood Method
#'
#' Computation of the log-likelihood of an extreme value mixture model (thus also \code{AIC} and \code{BIC} are available).
#'
#' @param object an object of class \code{evmm}.
#' @param ... additional parameters for compatibility.
#'
#' @importFrom stats logLik
#'
#' @return The log-likelihood of a model estimated with \code{extrememix}
#'
#' @examples logLik(rainfall_ggpd)
#' @method logLik evmm
#'
#' @export
#'
#' @rdname logLik
#'
logLik.evmm <- function(object,...){
  x <- object
  params <- apply(x$chain,2, median)
  if(sum(class(x) == "ggpd") == 1){ll <- sum(dggpd(x$data,params[1],params[2],params[3],params[4],params[5],log = T))
  attr(ll, "df") <- 5}
  if(sum(class(x) == "mgpd") == 1){k <- (ncol(x$chain)-3)/3
  ll <- sum(dmgpd(x$data,params[1],params[2],params[3],params[4:(4+k-1)],params[(4+k):(4+2*k-1)],params[(4+2*k):ncol(x$chain)],log = T))
  attr(ll, "df") <- (ncol(object$chain)-1)
  }
  attr(ll,"nobs") <- length(object$data)
  class(ll) <- "logLik"
  return(ll)
}
