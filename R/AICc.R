#' AICc model selection criterion
#'
#' Computation of the AICc for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name AICc
#' @return The AICc of a model estimated with \code{extrememix}
#' @export
AICc <- function (x, ...) {
  UseMethod("AICc", x)
}



#' @method AICc ggpd
#'@export
#' @rdname AICc
#'
AICc.ggpd <- function(x,...){
  AIC(x) + 60/(length(x$data)-6)}


#' @method AICc mgpd
#'@export
#' @rdname AICc
#'
AICc.mgpd <- function(x,...){
 AIC(x) + (2*ncol(x$chain)^2 + 2*ncol(x$chain))/(length(x$data)- ncol(x$chain)-1)
}