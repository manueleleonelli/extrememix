#' Quantiles Method
#'
#' Compute and plot quantiles
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name quantile
#' @export
quantile <- function (x, ...) {
  UseMethod("quantile", x)
}



#' @method quantile ggpd
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#'@export
#' @rdname quantile
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
quantile.ggpd <- function(x,values = NULL, cred = 0.95, ...){
   if(is.null(values)) {values <- c(0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.9925,0.995,0.9975,0.999,0.99925,0.9995)}
   out <- c_quant_ggpd(x$chain,values)
   mean <- round(apply(out,2,mean),2)
   lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
   upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
   empirical <- unname(stats::quantile(x$data,values))
   quantiles <- cbind(values,mean, lower, upper,empirical)
   colnames(quantiles) <- c("quantiles","estimate","lower_ci","upper_ci","empirical")
   output <- list(quantiles = quantiles, data = x$data, complete = out)
   class(output) <- "quantile"
   return(output)
  }


#' @method quantile mgpd
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#'@export
#' @rdname quantile
#'
#'@param values vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
quantile.mgpd <- function(x,values = NULL, cred = 0.95, ...){
   if(is.null(values)) {values <- c(0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.9925,0.995,0.9975,0.999,0.99925,0.9995)}
   k <- (ncol(x$chain)-3)/3;
   gpd <- x$chain[,1:3]
   mu <- x$chain[,4:(4+k-1)]
   eta <- x$chain[,(4+k):(4+2*k-1)]
   w <- x$chain[,(4+2*k):ncol(x$chain)]
   out <- c_quant_mgpd(gpd,mu,eta,w,values)
   mean <- round(apply(out,2,mean),2)
   lower <- round(apply(out,2,function(x)sort(x)[round(((1-cred)/2)*nrow(out))]),2)
   upper <- round(apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))]),2)
   empirical <- unname(stats::quantile(x$data,values))
   quantiles <- cbind(values,mean, lower, upper,empirical)
   colnames(quantiles) <- c("quantiles","estimate","lower_ci","upper_ci","empirical")
   output <- list(quantiles = quantiles, data = x$data, complete = out)
   class(output) <- "quantile"
   return(output)
}


