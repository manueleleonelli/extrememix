#' Log-likelihood method
#' @method logLik ggpd
#'@export
#' @rdname logLik
#'
logLik.ggpd <- function(object,...){
  x <- object
  params <- apply(x$chain,2, median)
  ll <- sum(dggpd(x$data,params[1],params[2],params[3],params[4],params[5],log = T))
  attr(ll, "df") <- 5
  attr(ll,"nobs") <- length(object$data)
  class(ll) <- "logLik"
  return(ll)
}


#' @method logLik mgpd
#'@export
#' @rdname logLik
#'
logLik.mgpd <- function(object,...){
  x <- object
  params <- apply(x$chain,2, mean)
  k <- (ncol(x$chain)-3)/3
  ll <- sum(dmgpd(x$data,params[1],params[2],params[3],params[4:(4+k-1)],params[(4+k):(4+2*k-1)],params[(4+2*k):ncol(x$chain)],log = T))
  attr(ll, "df") <- (ncol(object$chain)-1)
  attr(ll,"nobs") <- length(object$data)
  class(ll) <- "logLik"
  return(ll)
  }
