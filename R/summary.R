
#' Printing Methods
#'
#' Collection of printing methods for various objects created by \code{extrememix}.
#'
#' @param x an object created by \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name print
#' @returns A printed output of a model estimated with \code{extrememix}.
NULL


#' @rdname print
#'
#' @export
print.evmm <- function(x, ...) {
  if(sum(class(x) == "ggpd") == 1) cat("EVMM with Gamma bulk. LogLik", logLik(x), "\nxi estimated as ", median(x$chain[,1]), "\nProbability of unbounded distribution ", upper_bound(x)$prob)
  if(sum(class(x) == "mgpd") == 1) cat("EVMM with", (ncol(x$chain)-3)/3, "Mixtures of Gamma bulk. LogLik", logLik(x), "\nxi estimated as ", median(x$chain[,1]), "\nProbability of unbounded distribution ", upper_bound(x)$prob)
  invisible(x)
}


#' Summary Method
#'
#' Posterior estimates and credibility intervals for the parameters of extreme value mixture models.
#'
#' @param object an object of class \code{evmm}.
#' @param ... additional parameters (compatibility).
#'
#' @returns A printed summary of a model estimated with \code{extrememix} or any quantity associated with it.
#' @export

summary.evmm <- function(object,...){
  x <- object
  mean <- round(apply(x$chain, 2, median),2)
  upper <- round(apply(x$chain,2, function(x) sort(x)[round(0.975*length(x))]),2)
  lower <- round(apply(x$chain,2, function(x) sort(x)[round(0.025*length(x))]),2)
  if(sum(class(x) == "ggpd") == 1){
    names <- c("xi", "sigma","u","mu","eta")
    out <- data.frame(mean,lower,upper)
    rownames(out) <- names
    colnames(out) <- c("estimate","lower_ci","upper_ci")
    return(out)
  }
  if(sum(class(x) == "mgpd") == 1){
    k <- (ncol(x$chain) -3)/3
    if(k ==2){mu <- c("mu1","mu2"); eta <- c("eta1","eta2"); w <- c("w1","w2")}
    if(k ==3){mu <- c("mu1","mu2","mu3"); eta <- c("eta1","eta2","eta3"); w <- c("w1","w2","w3")}
    if(k ==4){mu <- c("mu1","mu2","mu3","mu4"); eta <- c("eta1","eta2","eta3","eta4"); w <- c("w1","w2","w3","w4")}
    names <- c("xi", "sigma","u",mu,eta,w)
    out <- data.frame(mean,lower,upper)
    rownames(out) <- names
    colnames(out) <- c("estimate","lower_ci","upper_ci")
    return(out)
  }
}

#'@rdname print
#'@export
print.summary.ggpd <- function(x,...){
  x
}

#'@rdname print
#'@export
print.quantile <- function(x, ...) {
  print(x$quantiles)
}

#'@rdname print
#'@export
print.return_level <- function(x, ...) {
  print(x$quantiles)
}


#'@rdname print
#'@export
print.VaR <- function(x, ...) {
  print(x$quantiles)
}

#'@rdname print
#'@export
print.ES <- function(x, ...) {
  print(x$quantiles)
}


#'@rdname print
#'@export
print.upper_bound <- function(x, ...) {
  if(x$prob < 1){
    cat("Probability of unbounded distribution: ", x$prob, "\nEstimated upper bound at ", round(median(x$bound),2), " with probability ", 1-x$prob, "\n Credibility interval at ", x$cred, "%: (", round(sort(x$bound)[round((1-x$cred)/2*length(x$bound))],2),",", round(sort(x$bound)[round((x$cred + (1-x$cred)/2)*length(x$bound))],2) ,")" )
    invisible(x)
  }
  else{
    cat("Probability of unbounded distribution: 1")
  }
}


