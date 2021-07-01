#' Predictive Distribution
#'
#' Plot of the predictive distribution of an extreme value mixture model
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name pred
#' @export
pred <- function (x, ...) {
  UseMethod("pred", x)
}



#' @method pred ggpd
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#'@export
#' @rdname pred
#'
#'@param x_axis vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
pred.ggpd <- function(x,x_axis = seq(0,max(x$data),0.01), cred = 0.95, ...){
  ..density.. <- NULL
  out <- c_pred_ggpd(x_axis,x$chain)
  mean <- apply(out,2, mean)
  lower <- apply(out,2, function(x) sort(x)[round(((1-cred)/2)*nrow(out))])
  upper <- apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))])
  data <- data.frame(x_axis,mean,lower,upper)
  ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.4) + theme_bw() + xlab("x") + ylab("Predictive Distribution") 
}



#' @method pred mgpd
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#'@export
#' @rdname pred
#'
#'@param x_axis vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
pred.mgpd <- function(x,x_axis = seq(0,max(x$data),0.01), cred = 0.95, ...){
  ..density.. <- NULL
  k <- (ncol(x$chain)-3)/3
  out <- c_pred_mgpd(x_axis,x$chain[,1],x$chain[,2],x$chain[,3],x$chain[,4:(4+k-1)], x$chain[, (4+k):(4+2*k-1)], x$chain[,(4+2*k):ncol(x$chain)])
  mean <- apply(out,2, mean)
  lower <- apply(out,2, function(x) sort(x)[round(((1-cred)/2)*nrow(out))])
  upper <- apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))])
  data <- data.frame(x_axis,mean,lower,upper)
  ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.4) + theme_bw() + xlab("x") + ylab("Predictive Distribution") 
}





