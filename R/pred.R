#' Predictive Distribution
#'
#' Plot of the predictive distribution of an extreme value mixture model.
#'
#' Consider an extreme value mixture model \eqn{f(y|\theta)} and suppose a sample \eqn{(\theta^{(1)},\dots,\theta^{(S)})} from the posterior distribution is available. The predictive distribution at the point \eqn{y} is estimated as
#' \deqn{\frac{1}{S}\sum_{s=1}^Sf(y|\theta^{(s)})}
#'
#' @references do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
#'
#' @param x the output of a model estimated with \code{extrememix}.
#' @param ... additional arguments for compatibility.
#' @name pred
#'
#' @return A plot of the estimate of the predictive distribution together with the data histogram.
#'
#' @examples pred(rainfall_ggpd)
#'
#' @export
pred <- function (x, ...) {
  UseMethod("pred", x)
}



#' @method pred evmm
#' @importFrom ggplot2 ggplot geom_line geom_ribbon theme_bw
#' @import ggplot2
#' @importFrom stats density
#'@export
#' @rdname pred
#'
#'@param x_axis vector of points where to estimate the predictive distribution.
#'@param cred amplitude of the posterior credibility interval.
#'@param xlim limits of the x-axis.
#'@param ylim limits of the y-axis.
pred.evmm <- function(x, x_axis = seq(min(x$data), max(x$data), length.out = 1000), cred = 0.95, xlim = c(min(x$data), max(x$data)), ylim = NULL, ...) {
  if(sum(class(x) == "ggpd") == 1) {
    out <- c_pred_ggpd(x_axis, x$chain)
  }
  if(sum(class(x) == "mgpd") == 1) {
    k <- (ncol(x$chain) - 3) / 3
    out <- c_pred_mgpd(x_axis, x$chain[,1], x$chain[,2], x$chain[,3], x$chain[,4:(4+k-1)], x$chain[, (4+k):(4+2*k-1)], x$chain[,(4+2*k):ncol(x$chain)])
  }

  mean <- apply(out, 2, median)
  lower <- apply(out, 2, function(x) sort(x)[round(((1 - cred) / 2) * nrow(out))])
  upper <- apply(out, 2, function(x) sort(x)[round((cred + (1 - cred) / 2) * nrow(out))])

  data <- data.frame(x_axis, mean, lower, upper)

  base_plot <- ggplot(data, aes(x = x_axis, y = mean)) +
    geom_histogram(data = data.frame(x = x$data), aes(x, y = after_stat(density)),
                   binwidth = 2 * length(x$data)^(-1/3) * (sort(x$data)[round(0.75 * length(x$data))] - sort(x$data)[round(0.25 * length(x$data))]),
                   colour = "black", fill = "#a6cee3", alpha = 0.6) +  # Blue fill for histogram
    geom_line(size = 1.2, color = "#e75480") +  # Red-pink line for predictive density
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#ff9999", alpha = 0.3) +  # Red-pink ribbon for credibility interval
    theme_minimal(base_size = 15) +  # Use minimal theme for clean look
    theme(
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.grid.major = element_line(color = "gray80", size = 0.5),  # Lighten major gridlines
      axis.text = element_text(color = "gray20"),  # Darken axis text for better readability
      axis.title = element_text(size = 14, face = "bold"),  # Make axis titles larger and bold
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  # Center and enhance title
    ) +
    xlab("x") +
    ylab("Predictive Distribution") +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)  # Allow zooming without distorting the axis

  if(is.null(ylim)) {
    return(base_plot)
  } else {
    return(base_plot + coord_cartesian(ylim = ylim, xlim = xlim, expand = FALSE))
  }
}
