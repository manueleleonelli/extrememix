#' Return Levels
#'
#' Computation of the return levels  for an extreme value mixture model
#'
#' A return level at \eqn{T} units of time is defined as the \eqn{1-1/T} quantile.
#'
#' @param x the output of a model estimated with \code{extrememix}
#' @param ... additional arguments for compatibility.
#' @seealso  \code{\link{ES}}, \code{\link{quant}}, \code{\link{VaR}}
#'
#' @references do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
#'
#'
#' @return A list with the following entries: \itemize{
#' \item \code{quantiles}: a matrix containing the estimated return levels, the posterior credibility intervals and the empirical estimate.
#' \item \code{data}: the dataset used to estimate the return levels.
#' \item \code{complete}: a matrix with the return levels for each value in the posterior sample.
#' }
#'
#' @examples return_level(rainfall_ggpd)
#'
#' @name return_level
#' @export
return_level<- function (x, ...) {
  UseMethod("return_level", x)
}


#' @method return_level evmm
#'@export
#' @rdname return_level
#'
#'@param values numeric vector of values of which to compute the value at risk.
#'@param cred amplitude of the posterior credibility interval.
return_level.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(20,25,30,40,50,60,70,80,90,100,150,200,250)}
  quantil <- 1-1/values
  out <- quant(x,quantil)
  out$quantiles[,1] <- values
  colnames(out$quantiles)[1] <- "Level"
  class(out) <- "return_level"
  return(out)
}
