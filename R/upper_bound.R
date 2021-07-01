#' Upper Bound Method
#'
#' Compute and plot return levels
#'
#' @param x the output of a model estimated with \code{extremix}
#' @param ... for compatibility
#' @name upper_bound
#' @export
 upper_bound <- function (x, ...) {
  UseMethod("upper_bound", x)
}


#' @method upper_bound evmm
#'@export
#' @rdname upper_bound
#'
#'@param cred amplitude of the posterior credibility interval
upper_bound.evmm <- function(x, cred = 0.95, ...){
  data <- x$chain[x$chain[,1] < 0,]
  subs <- apply(data,1,function(x) x[3] - x[2]/(x[1]*(1+x[1])))
  out <- list(bound = subs, prob = 1- nrow(data)/nrow(x$chain), cred = cred)
  class(out) <- "upper_bound"
  return(out)
}
