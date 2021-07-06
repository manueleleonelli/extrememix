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
#'@param x ciao
#'@param x_axis vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
#'@param xlim ciao
#'@param ylim ciao
pred.ggpd <- function(x,x_axis = seq(min(x$data),max(x$data),length.out = 1000), cred = 0.95, xlim = c(min(x$data),max(x$data)), ylim = NULL, ...){
  ..density.. <- NULL
  out <- c_pred_ggpd(x_axis,x$chain)
  mean <- apply(out,2, median)
  lower <- apply(out,2, function(x) sort(x)[round(((1-cred)/2)*nrow(out))])
  upper <- apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))])
  data <- data.frame(x_axis,mean,lower,upper)
  if(is.null(ylim)){ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.2) + theme_bw() + xlab("x") + ylab("Predictive Distribution") + coord_cartesian(xlim = xlim, expand = FALSE)}
  else{ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.2) + theme_bw() + xlab("x") + ylab("Predictive Distribution") +coord_cartesian(ylim = ylim, xlim= xlim, expand = FALSE) }

}



#' @method pred mgpd
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#'@export
#' @rdname pred
#'
#'@param x ciao
#'@param x_axis vector of points where to estimate the predictive distribution
#'@param cred amplitude of the posterior credibility interval
#'@param xlim ciao
#'@param ylim ciao
pred.mgpd <- function(x,x_axis = seq(min(x$data),max(x$data),length.out = 1000), cred = 0.95, xlim = c(min(x$data),max(x$data)), ylim = NULL, ...){
  ..density.. <- NULL
  k <- (ncol(x$chain)-3)/3
  out <- c_pred_mgpd(x_axis,x$chain[,1],x$chain[,2],x$chain[,3],x$chain[,4:(4+k-1)], x$chain[, (4+k):(4+2*k-1)], x$chain[,(4+2*k):ncol(x$chain)])
  mean <- apply(out,2, median)
  lower <- apply(out,2, function(x) sort(x)[round(((1-cred)/2)*nrow(out))])
  upper <- apply(out,2, function(x) sort(x)[round((cred+(1-cred)/2)*nrow(out))])
  data <- data.frame(x_axis,mean,lower,upper)
  if(is.null(ylim))ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.2) + theme_bw() + xlab("x") + ylab("Predictive Distribution") +coord_cartesian(xlim= xlim, expand = FALSE)
  else  ggplot(data,aes(x_axis,mean)) + geom_histogram(data = data.frame(x=x$data), aes(x, y = ..density..), binwidth = 2*length(x$data)^(-1/3)*(sort(x$data)[round(0.75*length(x$data))]-sort(x$data)[round(0.25*length(x$data))]), colour="black", fill="white") + geom_line() + geom_ribbon(aes(ymin=lower,ymax = upper),alpha =0.2) + theme_bw() + xlab("x") + ylab("Predictive Distribution") +coord_cartesian(ylim = ylim, xlim= xlim, expand = FALSE)
  }





