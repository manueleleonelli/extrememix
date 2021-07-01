#' DIC model selection criterion
#'
#' Computation of the WAIC for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name WAIC
#' @return The WAIC of a model estimated with \code{extrememix}
#' @export
WAIC <- function (x, ...) {
  UseMethod("WAIC", x)
}



#' @method WAIC ggpd
#'@export
#' @rdname WAIC
#'
WAIC.ggpd <- function(x,...){
  WAIC_ggpd(x$chain,x$data)
}



