#' BIC model selection criterion
#'
#' Computation of the BIC for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name BIC
#' @return The BIC of a model estimated with \code{extrememix}
#' @export
BIC <- function (x, ...) {
  UseMethod("BIC", x)
}



#' @method BIC ggpd
#'@export
#' @rdname BIC
#'
BIC.ggpd <- function(x,...){
 params <- apply(x$chain,2, mean)
 -2*sum(dggpd(x$data,params[1],params[2],params[3],params[4],params[5],log = T)) +5*log(length(x$data))
}


#' @method BIC mgpd
#'@export
#' @rdname BIC
#'
BIC.mgpd <- function(x,...){
  params <- apply(x$chain,2, mean)
  -2*logLik(x) +(ncol(x$chain)-1)*log(length(x$data))
}
