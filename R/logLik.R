#' Log-likelihood method
#'
#' Computation of the log-likelihood for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name logLik
#' @return The log-likelihood of a model estimated with \code{extrememix}
#' @export
logLik <- function (x, ...) {
  UseMethod("logLik", x)
}



#' @method logLik ggpd
#'@export
#' @rdname logLik
#'
logLik.ggpd <- function(x,...){
  params <- apply(x$chain,2, mean)
  sum(dggpd(x$data,params[1],params[2],params[3],params[4],params[5],log = T))
}


#' @method logLik mgpd
#'@export
#' @rdname logLik
#'
logLik.mgpd <- function(x,...){
  params <- apply(x$chain,2, mean)
  k <- (ncol(x$chain)-3)/3
  sum(dmgpd(x$data,params[1],params[2],params[3],params[4:(4+k-1)],params[(4+k):(4+2*k-1)],params[(4+2*k):ncol(x$chain)],log = T))
}
