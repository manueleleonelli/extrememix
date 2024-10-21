
#' Plot Methods for Summaries
#'
#' Plotting methods for objects created with \code{quant}, \code{ES}, \code{return_level} or \code{VaR}.
#'
#' Two types of plot can be output: either a line plot in the case the functions \code{quant}, \code{ES}, \code{return_level} or \code{VaR} where called with more than one value for the input \code{values}, or an histogram otherwise.
#'
#' @param x an object of class \code{quant}, \code{ES}, \code{return_level} or \code{VaR}.
#' @param ylim limits of the y-axis.
#' @param ... additional parameters  for compatibility.
#' @returns Appropriate plots for quantities computed with \code{extrememix}.
#'
#' @name plot_summaries
#'
#' @examples
#' plot(return_level(rainfall_ggpd)) ## for line plot
#' plot(return_level(rainfall_ggpd, values = 100)) ## for histogram
#'
#'
NULL



#' @rdname plot_summaries
#' @export
plot.quant <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    quantiles <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(quantiles,estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Probability",  y= "Quantile", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else {ggplot(data.frame(x$quantiles),aes(quantiles,estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Probability",  y= "Quantile", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) + coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
     x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("quantile")
  }
}

#' @rdname plot_summaries
#' @export
plot.return_level <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Levels",  y= "Return Levels", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Levels",  y= "Return Levels", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)+coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("Return Level")
  }
}

#' @rdname plot_summaries
#' @export
plot.VaR <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    VaR_Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=VaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Var_Levels",  y= "VaR", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=VaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Var_Levels",  y= "VaR", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) + coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("VaR")
  }
}

#' @rdname plot_summaries
#' @export
plot.ES <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    ES_Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=ES_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "ES_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=ES_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "ES_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) +coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("ES")
  }
}



#' Plot Upper Bounds
#'
#' Plotting method for the posterior distribution of the upper bound. No plot is reported if the posterior sample of xi has only positive values (unbounded distribution).
#'
#' @param x an object of class \code{upper_bound}.
#' @param xlim limits of the x-axis.
#' @param ... additional parameters for compatibility.
#'
#' @examples plot(upper_bound(rainfall_ggpd))
#'
#'@returns A histogram for the posterior estimated upper bound of the distribution.
#' @export

plot.upper_bound <- function(x, xlim = c(min(x$bound),max(x$bound)),...) {
  IQR <-  NULL
  if (x$prob == 1){stop("No plot created since probability of unbounded distribution is 1")}
  else{
    x.bound <- x$bound[x$bound >= xlim[1] & x$bound <= xlim[2]]
    ggplot(data.frame(x.bound), aes(x=x.bound)) +    geom_histogram(aes(y=after_stat(density)), binwidth =  2*length(x.bound)^(-1/3)*IQR(x.bound), colour="black", fill="white") + theme_bw()+ xlab(cat("Upper Bound, with probability ", 1 - x$prob))
  }
}

#' Plot of Extreme Value Mixture Models
#'
#' Plotting method for objects of class \code{evmm} giving an overview of an estimated model.
#'
#' The \code{plot} method for objects of class \code{evmm} reports four plots: \itemize{
#' \item An histogram of the posterior distribution of xi.
#' \item An histogram of the posterior distribution of sigma.
#' \item A line plot of the estimated quantiles, their posterior credibility interval, and the empirical ones.
#' \item A plot of the predictive distribution together with the data histogram.
#' }
#'
#' @examples plot(rainfall_ggpd)
#'
#' @returns Plots of a model estimated with \code{extrememix}.
#' @import gridExtra
#' @param x an object of class \code{evmm}.
#' @param ... additional parameters for compatibility.
#'
#' @export

plot.evmm <- function(x, ...) {
  X1 <-X2 <- IQR <- NULL
  data <- data.frame(x$chain)
  colnames(data)[1] <- "X1"
  colnames(data)[2] <- "X2"
  p1 <- ggplot(data, aes(x = X1)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*nrow(x$chain)^(-1/3)*IQR(x$chain[,1]), colour="black", fill="white") + labs(x = "xi") + theme_bw() +geom_vline(xintercept=median(x$chain[,1]), linetype="dashed", color = "red")
  p2 <- ggplot(data, aes(x = X2)) + geom_histogram(aes(y=after_stat(density)), binwidth = 2*length(x$chain)^(-1/3)*IQR(x$chain[,2]), colour="black", fill="white") + labs(x = "sigma") + theme_bw() +geom_vline(xintercept=median(x$chain[,2]), linetype="dashed",  color = "red")
  p3 <- plot(quant(x),...)
  p4 <- pred(x,...)
  grid.arrange(p1,p2,p3,p4,ncol = 2, nrow =2)
}
