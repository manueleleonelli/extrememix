#' Value-at-Risk
#'
#' Computation of the Value-at-Risk for an extreme value mixture model.
#'
#' The Value-at-Risk for level q\% is defined as the 1-q/100 quantile of the model.
#'
#' @param x the output of a model estimated with \code{extrememix}
#' @param ... additional arguments for compatibility.
#'
#' @seealso  \code{\link{ES}}, \code{\link{quant}}, \code{\link{return_level}}
#'
#' @references Lattanzi, Chiara, and Manuele Leonelli. "A changepoint approach for the identification of financial extreme regimes." Brazilian Journal of Probability and Statistics.
#'
#'
#' @return A list with the following entries: \itemize{
#' \item \code{quantiles}: a matrix containing the estimated value at risk, the posterior credibility intervals and the empirical estimate.
#' \item \code{data}: the dataset used to estimate the value at risk.
#' \item \code{complete}: a matrix with the value at risk for each value in the posterior sample.
#' }
#'
#' @examples VaR(rainfall_ggpd)

#' @name VaR
#' @export
VaR <- function (x, ...) {
  UseMethod("VaR", x)
}


#' @method VaR evmm
#'@export
#' @rdname VaR
#'
#'@param values numeric vector of values of which to compute the value at risk.
#'@param cred amplitude of the posterior credibility interval.
VaR.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quantil <- 1-values/100
  out <- quant(x,quantil)
  out$quantiles[,1] <- values
  colnames(out$quantiles)[1] <- "VaR_Level"
  class(out) <- "VaR"
  return(out)
}
