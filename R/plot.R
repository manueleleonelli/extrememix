#' Plot quantiles
#'
#' @param x an object of class \code{quant}.
#' @param ylim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.quant <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
  quantiles <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
  colors <- c("est" = "black", "emp" = "red")
  if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(quantiles,estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Probability",  y= "Quantile", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
  else {ggplot(data.frame(x$quantiles),aes(quantiles,estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Probability",  y= "Quantile", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) + coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    ..density.. <- x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("quantile")
  }
}

#' Plot return levels
#'
#' @param x an object of class \code{return_level}.
#' @param ylim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.return_level <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Levels",  y= "Return Levels", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Levels",  y= "Return Levels", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)+coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    ..density.. <- x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("Return Level")
  }
}

#' Plot VaRs
#'
#' @param x an object of class \code{VaR}.
#' @param ylim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.VaR <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    VaR_Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=VaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Var_Levels",  y= "VaR", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=VaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "Var_Levels",  y= "VaR", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) + coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    ..density.. <- x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("VaR")
  }
}

#' Plot ESs
#'
#' @param x an object of class \code{ES}.
#' @param ylim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.ES <- function(x, ylim = NULL, ...) {
  if(nrow(x$quantiles) > 1){
    ES_Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=ES_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "ES_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=ES_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "ES_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) +coord_cartesian(ylim = ylim, expand = FALSE)}
  }
  else{
    ..density.. <- x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$complete)^(-1/3)*(sort(x$complete)[round(0.75*length(x$complete))]-sort(x$complete)[round(0.25*length(x$complete))]), colour="black", fill="white") +theme_bw()+xlab("ES")
  }
}


#' Plot TVaRs
#'
#' @param x an object of class \code{TVaR}.
#' @param ylim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.TVaR <- function(x, ylim = NULL, ...) {
  IQR <- NULL
  if(nrow(x$quantiles) > 1){
    TVaR_Level <- estimate <- lower_ci <- upper_ci <- empirical <- NULL
    colors <- c("est" = "black", "emp" = "red")
    if(is.null(ylim)){ggplot(data.frame(x$quantiles),aes(x=TVaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "TVaR_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors)}
    else{ggplot(data.frame(x$quantiles),aes(x=TVaR_Level,y=estimate))+ geom_line(aes(color = "est")) + geom_ribbon(aes(ymin=lower_ci,ymax = upper_ci),alpha =0.2) + theme_bw() + labs(x= "TVaR_Levels",  y= "ES", color = "") + geom_line(aes(y=empirical, col='emp'))+ scale_color_manual(values = colors) +coord_cartesian(ylim = ylim, expand = FALSE) }
  }
  else{
    ..density.. <- x.complete <-  NULL
    ggplot(data.frame(x$complete), aes(x = x.complete)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$complete)^(-1/3)*IQR(x$complete), colour="black", fill="white") +theme_bw()+xlab("TVaR")
  }
}

#' Plot Upper Bounds
#'
#' @param x an object of class \code{upper_bound}.
#' @param xlim ciao
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.upper_bound <- function(x, xlim = c(min(x$bound),max(x$bound)),...) {
  IQR <- ..density.. <- NULL
if (x$prob == 1){stop("No plot created since probability of unbounded distribution is 1")}
  else{
    x.bound <- x$bound[x$bound >= xlim[1] & x$bound <= xlim[2]]
    ggplot(data.frame(x.bound), aes(x=x.bound)) +    geom_histogram(aes(y=..density..), binwidth =  2*length(x.bound)^(-1/3)*IQR(x.bound), colour="black", fill="white") + theme_bw()+ xlab(cat("Upper Bound, with probability ", 1 - x$prob)) 
  }
}

#' Plot quantiles
#' @import gridExtra
#' @param x an object of class \code{evmm}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

plot.evmm <- function(x, ...) {
  X1 <- ..density.. <- X2 <- IQR <- NULL
  p1 <- ggplot(data.frame(x$chain), aes(x = X1)) + geom_histogram(aes(y=..density..), binwidth = 2*nrow(x$chain)^(-1/3)*IQR(x$chain[,1]), colour="black", fill="white") + labs(x = "xi") + theme_bw() +geom_vline(xintercept=median(x$chain[,1]), linetype="dashed", 
                                                                                                                                         color = "red")
  p2 <- ggplot(data.frame(x$chain), aes(x = X2)) + geom_histogram(aes(y=..density..), binwidth = 2*length(x$chain)^(-1/3)*IQR(x$chain[,2]), colour="black", fill="white") + labs(x = "sigma") + theme_bw() +geom_vline(xintercept=median(x$chain[,2]), linetype="dashed", 
                                                                                                                                                                                                      color = "red")
  p3 <- plot(quant(x),...)
  p4 <- pred(x,...)
  grid.arrange(p1,p2,p3,p4,ncol = 2, nrow =2)
}
