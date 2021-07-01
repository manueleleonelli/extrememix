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


#' @method WAIC mgpd
#'@export
#' @rdname WAIC
#'
WAIC.mgpd <- function(x,...){
  k <- (ncol(x$chain)-3)/3
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  WAIC_mgpd(gpd,mu,eta,w,x$data)
}
  

