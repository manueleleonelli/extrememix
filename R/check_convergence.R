#' Convergence Assessment of MCMC Algorithms
#'
#' Plot of the traceplot and autocorrelation function for the 0.99 quantile from the posterior sample.
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name check_convergence
#' @export
#' @return Two plots to check if the estimation with \code{fggpd} and \code{mgpd} converged: traceplot and autocorrelation plot for the 99th quantile of the posterior density.
#' @examples check_convergence(rainfall_ggpd)
check_convergence <- function (x, ...) {
  UseMethod("check_convergence", x)
}

#' @method check_convergence evmm
#' @import ggplot2
#' @import gridExtra
#'@export
#' @rdname check_convergence
#'
check_convergence.evmm <- function(x, ...){
  quantil <- quant(x,0.99)$complete
  p1 <- ggplot(data.frame(quantil), aes(x = 1:length(quantil),y = quantil)) + geom_line() + labs(x ="Index", y = "Quantile") + theme_bw()
  p2 <- suppressWarnings(ggacf(quantil))
  grid.arrange(p1,p2,ncol=2)
}

#' @importFrom stats qnorm
ggacf <- function(series) {
  acf <- lag <- NULL
  significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(series)))
  a<-acf(series, plot=F)
  a.2<-with(a, data.frame(lag, acf))
  g<- ggplot(a.2[-1,], aes(x=lag,y=acf)) +
    geom_bar(stat = "identity", position = "identity") + xlab('Lag') + ylab('ACF') +
    geom_hline(yintercept=c(significance_level,-significance_level), lty=3);

  # fix scale for integer lags
  if (all(a.2$lag%%1 == 0)) {
    g<- g + scale_x_discrete(limits = seq(1, max(a.2$lag))) + theme_bw();
  }
  return(g);
}
