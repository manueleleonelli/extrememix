#' DIC model selection criterion
#'
#' Computation of the DIC for an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name DIC
#' @return The DIC of a model estimated with \code{extrememix}
#' @export
DIC <- function (x, ...) {
  UseMethod("DIC", x)
}



#' @method DIC ggpd
#'@export
#' @rdname DIC
#'
DIC.ggpd <- function(x,...){
 DIC_ggpd(x$chain,x$data)
}


#' @method DIC mgpd
#'@export
#' @rdname DIC
#'
DIC.mgpd <- function(x,...){
  k <- (ncol(x$chain)-3)/3
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  DIC_mgpd(gpd,mu,eta,w,x$data)
  
  
#  d <- -2*logLik(x)
#  k <- (ncol(x$chain)-3)/3
#  bard <- 0
#  for(i in 1:nrow(x$chain)){
#    bard <- bard -2* sum(dmgpd(x$data,x$chain[i,1],x$chain[i,2],x$chain[i,3],x$chain[i,4:(4+k-1)],x$chain[i,(4+k):(4+2*k-1)],x$chain[i,(4+2*k):ncol(x$chain)],log = T))
#  }
#  bard <- (1/nrow(x$chain))*bard
#  return(2*bard - d)
}



