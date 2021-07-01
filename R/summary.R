#' Print a ggpd
#'
#' @param x an object of class \code{ggpd}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.ggpd <- function(x, ...) {
  cat("EVMM with Gamma bulk. LogLik", logLik(x), "\nxi estimated as ", mean(x$chain[,1]), "\nProbability of unbounded distribution ", upper_bound(x)$prob)
  invisible(x)
}

#' Summary of a ggpd
#' @param object an object of class \code{ggpd}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

summary.ggpd <- function(object,...){
  x <- object
 mean <- round(apply(x$chain, 2, mean),2)
 upper <- round(apply(x$chain,2, function(x) sort(x)[round(0.975*length(x))]),2)
 lower <- round(apply(x$chain,2, function(x) sort(x)[round(0.025*length(x))]),2)
 names <- c("xi", "sigma","u","mu","eta")
 out <- data.frame(mean,lower,upper)
 rownames(out) <- names
 colnames(out) <- c("estimate","lower_ci","upper_ci")
 return(out)

}

#'@export
print.summary.ggpd <- function(x,...){
 x
}

#' Print a quantile object
#'
#' @param x an object of class \code{guantile}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.quantile <- function(x, ...) {
  print(x$quantiles)
}

#' Print a return_level object
#'
#' @param x an object of class \code{return_level}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.return_level <- function(x, ...) {
  print(x$quantiles)
}


#' Print a VaR object
#'
#' @param x an object of class \code{VaR}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.VaR <- function(x, ...) {
  print(x$quantiles)
}


#' Print a ES object
#'
#' @param x an object of class \code{ES}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.ES <- function(x, ...) {
  print(x$quantiles)
}

#' Print a TVaR object
#'
#' @param x an object of class \code{ES}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.TVaR <- function(x, ...) {
  print(x$quantiles)
}


#' Print a upper bound 
#'
#' @param x an object of class \code{upper_bound}.
#' @param ... additional parameters (compatibility).
#'
#' @details Info
#' @export

print.upper_bound <- function(x, ...) {
  if(x$prob < 1){
  cat("Probability of unbounded distribution: ", x$prob, "\nEstimated upper bound at ", round(mean(x$bound),2), " with probability ", 1-x$prob, "\n Credibility interval at ", x$cred, "%: (", round(sort(x$bound)[round((1-x$cred)/2*length(x$bound))],2),",", round(sort(x$bound)[round((x$cred + (1-x$cred)/2)*length(x$bound))],2) ,")" )
  invisible(x)
  }
  else{
    cat("Probability of unbounded distribution: 1")
  }
}


