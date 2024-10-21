#' Estimated Quantiles
#'
#' Computation of posterior quantiles for an extreme value mixture model
#'
#' For a random variable \eqn{X} the p-quantile is the value \eqn{x} such that \eqn{P(X>x)=1-p}. For an extreme value mixture model this can be computed as
#' \deqn{x = u +\frac{\sigma}{\xi}((1-p^*)^{-\xi}-1),}
#' where \deqn{p^* = \frac{p-F_\textnormal{bulk}(u|\theta)}{1-F_\textnormal{bulk}(u|\theta)},}
#'and \eqn{F_\textnormal{bulk}} is the distribution function of the bulk, parametrized by \eqn{\theta}.
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name quant
#'
#'@references do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
#'
#'
#' @return A list with the following entries: \itemize{
#' \item \code{quantiles}: a matrix containing the quantiles, the posterior credibility intervals and the empirical estimate.
#' \item \code{data}: the dataset used to estimate the quantiles.
#' \item \code{complete}: a matrix with the quantiles for each value in the posterior sample.
#' }
#'
#' @examples quant(rainfall_ggpd)
#'
#'
#' @export
quant <- function (x, ...) {
  UseMethod("quant", x)
}



#' @method quant evmm
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#' @importFrom stats median
#'@export
#' @rdname quant
#'
#'@param values numeric vector of values of which to compute the quantile.
#'@param cred amplitude of the posterior credibility interval.
quant.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.9925,0.995)}
  if(sum(class(x) == "ggpd") == 1) {out <- c_quant_ggpd(x$chain,values)}
  if(sum(class(x) == "mgpd") == 1) {k <- (ncol(x$chain)-3)/3;
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  out <- c_quant_mgpd(gpd,mu,eta,w,values)}
  mean <- round(apply(out,2,median),2)
  lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
  upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
  empirical <- round(unname(stats::quantile(x$data,values)),2)
  quantiles <- cbind(values,mean, lower, upper,empirical)
  colnames(quantiles) <- c("quantiles","estimate","lower_ci","upper_ci","empirical")
  output <- list(quantiles = quantiles, data = x$data, complete = out)
  class(output) <- "quant"
  return(output)
}
