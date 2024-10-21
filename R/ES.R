#' Expected Shortfall
#'
#' Computation of the expected shortfall for an extreme value mixture model
#'
#' The expected shortfall is the expectation of a random variable conditional of being larger of a specific Value-at-Risk (quantile). For an extreme value mixture model this is equal to: \deqn{ES_p = \frac{VaR_p}{1-\xi} +\frac{\sigma-\xi u }{1-\xi}}
#'
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#'
#' @seealso  \code{\link{quant}}, \code{\link{return_level}}, \code{\link{VaR}}
#'
#' @references Lattanzi, Chiara, and Manuele Leonelli. "A changepoint approach for the identification of financial extreme regimes." Brazilian Journal of Probability and Statistics.
#'
#' @return A list with the following entries: \itemize{
#' \item \code{quantiles}: a matrix containing the estimated shortfall, the posterior credibility intervals and the empirical estimate.
#' \item \code{data}: the dataset used to estimate the expected shortfall.
#' \item \code{complete}: a matrix with the expected shortfall for each value in the posterior sample.
#' }
#'
#' @examples ES(rainfall_ggpd)
#'
#' @name ES
#' @export
ES <- function (x, ...) {
  UseMethod("ES", x)
}

median <- function(x){stats::quantile(x,0.5,type=2)}

#' @method ES evmm
#'@export
#' @rdname ES
#'
#'@param values numeric vector of values of which to compute the expected shortfall.
#'@param cred amplitude of the posterior credibility interval.
ES.evmm <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  if(sum(class(x) == "ggpd") == 1) {out <- c_es_ggpd(x$chain,quant)}
  if(sum(class(x) == "mgpd") == 1) {k <- (ncol(x$chain)-3)/3;
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  out <- c_es_mgpd(gpd,mu,eta,w,quant)}
  mean <- round(apply(out,2,median),2)
  lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
  upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
  empirical <- round(unname(stats::quantile(x$data,quant,type = 2)),2)
  for (i in 1:length(empirical)) empirical[i] <- mean(x$data[x$data>= empirical[i]])
  quantiles <- cbind(values,mean, lower, upper,empirical)
  colnames(quantiles) <- c("ES_Level","estimate","lower_ci","upper_ci","empirical")
  output <- list(quantiles = quantiles, data = x$data, complete = out)
  class(output) <- "ES"
  return(output)
}
