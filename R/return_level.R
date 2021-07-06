#' Return Levels Method
#'
#' Compute and plot return levels
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name return_level
#' @export
return_level<- function (x, ...) {
  UseMethod("return_level", x)
}


#' @method return_level evmm
#'@export
#' @rdname return_level
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
return_level.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(20,25,30,40,50,60,70,80,90,100,150,200,250)}
  quantil <- 1-1/values
  out <- quant(x,quantil)
  out$quantiles[,1] <- values
  colnames(out$quantiles)[1] <- "Level"
  class(out) <- "return_level"
  return(out)
}
