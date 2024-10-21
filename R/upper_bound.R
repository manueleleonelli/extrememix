#' Upper Bound
#'
#' Computation of the upper bound of the distribution
#'
#' For an extreme value mixture model with a shape parameter \eqn{xi < 0} the distribution is right-bounded with upper limit equal to \eqn{u-\sigma/\xi}.
#'
#' @return \code{upper_bound} returns a list with entries: \itemize{
#' \item \code{bound}: a sample from the posterior distribution of the upper limit of the model, taken over the posterior values of xi which are negative.
#' \item \code{prob}: the posterior probability that the distribution is unbounded.
#' \item \code{cred}: the requested amplitude of the posterior credibility intervals.
#' }
#'
#' @examples upper_bound(rainfall_ggpd)
#'
#' @references Coles, Stuart, et al. An introduction to statistical modeling of extreme values. Vol. 208. London: Springer, 2001.
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name upper_bound
#' @export
upper_bound <- function (x, ...) {
  UseMethod("upper_bound", x)
}


#' @method upper_bound evmm
#'@export
#' @rdname upper_bound
#'
#'@param cred amplitude of the posterior credibility interval.
upper_bound.evmm <- function(x, cred = 0.95, ...){
  data <- x$chain[x$chain[,1] < 0,]
  subs <- apply(data,1,function(x) x[3] - x[2]/(x[1]))
  out <- list(bound = subs, prob = 1- nrow(data)/nrow(x$chain), cred = cred)
  class(out) <- "upper_bound"
  return(out)
}
