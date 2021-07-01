#' AIC model selection criterion
#'
#' Computation of the AIC for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name AIC
#' @return The AIC of a model estimated with \code{extrememix}
#' @export
AIC <- function (x, ...) {
  UseMethod("AIC", x)
}



#' @method AIC ggpd
#'@export
#' @rdname AIC
#'
AIC.ggpd <- function(x,...){
  params <- apply(x$chain,2, mean)
  -2*sum(dggpd(x$data,params[1],params[2],params[3],params[4],params[5],log = T)) +5*2
}


#' @method AIC mgpd
#'@export
#' @rdname AIC
#'
AIC.mgpd <- function(x,...){
  params <- apply(x$chain,2, mean)
  -2*logLik(x) + (ncol(x$chain)-1)*2
}
