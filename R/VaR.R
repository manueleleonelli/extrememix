#' VaR Method
#'
#' Compute and plot VaRs
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name VaR
#' @export
VaR <- function (x, ...) {
  UseMethod("VaR", x)
}


#' @method VaR evmm
#'@export
#' @rdname VaR
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
VaR.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  out <- quantile(x,quant)
  out$quantiles[,1] <- values
  colnames(out$quantiles)[1] <- "VaR_Level"
  class(out) <- "VaR"
  return(out)
}
