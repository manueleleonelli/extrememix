#' Expected Shortfall Method
#'
#' Compute and plot Expected Shortfalls
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name ES
#' @export
ES <- function (x, ...) {
  UseMethod("ES", x)
}


#' @method ES ggpd
#'@export
#' @rdname ES
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
ES.ggpd <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <- c(0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  out <- c_es_ggpd(x$chain,quant)
  mean <- round(apply(out,2,mean),2)
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

#' @method ES mgpd
#'@export
#' @rdname ES
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
ES.mgpd <- function(x,values = NULL, cred = 0.95, ...){
  if(is.null(values)) {values <-c(0.5,0.75,1,1.5,2,2.5,3,4,5)}
  quant <- 1-values/100
  k <- (ncol(x$chain)-3)/3;
  gpd <- x$chain[,1:3]
  mu <- x$chain[,4:(4+k-1)]
  eta <- x$chain[,(4+k):(4+2*k-1)]
  w <- x$chain[,(4+2*k):ncol(x$chain)]
  out <- c_es_mgpd(gpd,mu,eta,w,quant)
  mean <- round(apply(out,2,mean),2)
  lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
  upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
  empirical <- round(unname(stats::quantile(x$data,quant, type = 2)))
  for (i in 1:length(empirical)) empirical[i] <- mean(x$data[x$data>= empirical[i]])
  quantiles <- cbind(values,mean, lower, upper,empirical)
  colnames(quantiles) <- c("ES_Level","estimate","lower_ci","upper_ci","empirical")
  output <- list(quantiles = quantiles, data = x$data, complete = out)
  class(output) <- "ES"
  return(output)
}
