#' VaR Method
#'
#' Compute and plot VaRs
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name check_convergence
#' @export
check_convergence <- function (x, ...) {
  UseMethod("check_convergence", x)
}

#' @method check_convergence evmm
#'@export
#' @rdname check_convergence
#'
check_convergence.evmm <- function(x, ...){
  quantil <- quant(x,0.99)$complete
  p1 <- ggplot(data.frame(quantil), aes(x = 1:length(quantil),y = quantil)) + geom_line() + labs(x ="Index", y = "Quantile") + theme_bw()
  p2 <- suppressWarnings(ggacf(quantil))
  grid.arrange(p1,p2,ncol=2)
}

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
