#' TVaR Method
#'
#' Compute and plot TVaRs
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name TVaR
#' @export
TVaR <- function (x, ...) {
  UseMethod("TVaR", x)
}


#' @method TVaR ggpd
#'@export
#' @rdname TVaR
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
TVaR.ggpd <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  out <- c_tvar_ggpd(x$chain,quant)
  mean <- round(apply(out,2,mean),2)
  lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
  upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
  empirical <- unname(stats::quantile(x$data,quant))
  for (i in 1:length(empirical)) empirical[i] <- mean(x$data[x$data>= empirical[i]]) + empirical[i]
  quantiles <- cbind(values,mean, lower, upper,empirical)
  colnames(quantiles) <- c("TVaR_Level","estimate","lower_ci","upper_ci","empirical")
  output <- list(quantiles = quantiles, data = x$data, complete = out)
  class(output) <- "TVaR"
  return(output)
}


#' @method TVaR mgpd
#'@export
#' @rdname TVaR
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
TVaR.mgpd <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  k <- (ncol(x$chain)-3)/3;
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  out <- c_tvar_mgpd(gpd,mu,eta,w,quant)
  mean <- round(apply(out,2,mean),2)
  lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
  upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
  empirical <- unname(stats::quantile(x$data,quant))
  for (i in 1:length(empirical)) empirical[i] <- mean(x$data[x$data>= empirical[i]]) + empirical[i]
  quantiles <- cbind(values,mean, lower, upper,empirical)
  colnames(quantiles) <- c("TVaR_Level","estimate","lower_ci","upper_ci","empirical")
  output <- list(quantiles = quantiles, data = x$data, complete = out)
  class(output) <- "TVaR"
  return(output)
}
