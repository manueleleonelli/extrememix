#' @import mixtools
#' @import evd
#' @importFrom stats quantile
#' @importFrom utils capture.output
#' @importFrom threshr ithresh
starting.values <- function(x, k){
  u <- unname(summary(ithresh(data = x,  u_vec =  stats::quantile(x, probs = seq(0.75, 0.95, by = 0.05)), n = 100)))[3]
  lower <- subset(x, x <= u)
  pargpd <- unname(fpot(x,stats::quantile(x,0.9),std.err = F)$estimate)
  invisible(capture.output(mix <- gammamixEM(x, k = k, epsilon = 1e-02)))
  mu <- unname(mix$gamma.pars[1,])*unname(mix$gamma.pars[2,])
  eta <- unname(mix$gamma.pars[1,])
  return(list(xi = max(-0.4,pargpd[2]),sigma = pargpd[1], u = u, mu = sort(mu) , eta = eta[order(mu)], w = rep(1/k,k)))
}


#' @importFrom stats uniroot
#' @importFrom stats qnorm
compute.var <- function(x,start,k){
  if(start$xi>=0 ){xi <- 0.025} else{xi <- 0.001}
  sigma <- uniroot(function(r) qnorm(0.01, start$sigma,r) - max(0.9*start$sigma,0) , c(0,max(x)))$root
  u <- uniroot(function(r) qnorm(0.01,start$u,r) - max(0.9*start$u,0), c(0,max(x)))$root
  mu <- rep(0,length(start$mu))
  for(i in 1:k){
    mu[i] <- uniroot(function(r) qnorm(0.01,start$mu[i],r) -  max(0.99*start$mu[i],0), c(0,max(x)))$root
  }
  w <- 0.1
  return(list(xi=xi, sigma = sigma, u = u, mu = mu, w = w))
}


check.input <- function(x, it, k, var, start, prior, thin, burn){
  if(any(x <=0)) stop("x must be positive")
  if(k%% 1 !=0 | k <= 0 ) stop("k must be a positive integer")
  if(k == 1) stop("use the function fggpd")
  if(k >4) stop("Maximum number of mixture components is 4")
  if(it %% 1 != 0 | it <= 0){stop('iterations must be a positive integer')}
  if(burn > it){stop('burn-in cannot be larger than iterations')}
  if(thin > (it-burn)){stop('thinning too large')}
  if(mode(var) != "list"){stop('var must be a list')}
  if(mode(start) != "list"){stop('start must be a list')}
  if(mode(prior) != "list"){stop('prior must be a list')}
  if(sum(start$w) - 1 > 0.000001) stop('weights must add to one')

  if(is.null(var[["xi"]])){stop('Variance for xi not provided')}
  if(is.null(var[["sigma"]])){stop('Variance for sigma not provided')}
  if(is.null(var[["u"]])){stop('Variance for u not provided')}
  if(is.null(var[["mu"]])){stop('Variance for mu not provided')}
  if(is.null(var[["w"]])){stop('Variance for w not provided')}

  if(is.null(start[["xi"]])){stop('Starting value for xi not provided')}
  if(is.null(start[["sigma"]])){stop('Starting value for sigma not provided')}
  if(is.null(start[["u"]])){stop('Starting value for u not provided')}
  if(is.null(start[["mu"]])){stop('Starting value for mu not provided')}
  if(is.null(start[["eta"]])){stop('Starting value for eta not provided')}
  if(is.null(start[["w"]])){stop('Starting value for w not provided')}

  if(is.null(prior[["u"]])){stop('Prior for u not provided')}
  if(is.null(prior[["mu_mu"]])){stop('Prior mean for mu not provided')}
  if(is.null(prior[["mu_eta"]])){stop('Prior shape for mu not provided')}
  if(is.null(prior[["eta_mu"]])){stop('Prior mean for eta not provided')}
  if(is.null(prior[["eta_eta"]])){stop('Prior shape for eta not provided')}
}


#' @importFrom stats uniroot
#' @importFrom stats qnorm
compute.prior <- function(x,start){
  eta_mu <- start$eta
  eta_eta <- rep(0.001, length(start$eta))
  prior_u <- c(start$u,  uniroot(function(r) qnorm(0.025,start$u,r) - stats::quantile(x,0.5), c(0,max(x)))$root)
  mu_mu <- start$mu
  mu_eta <- rep(0.001, length(start$eta))
  return(list(u = prior_u, mu_mu = mu_mu, mu_eta = mu_eta, eta_mu = eta_mu, eta_eta = eta_eta))
}


#' MGPD Estimation
#'
#' Fit of the MGPD model using an MCMC algorithm.
#'
#' Estimation of the MGPD is carried out using an adaptive block Metropolis-Hastings algorithm. As standard, the user needs to specify the data to use during estimation, the number of mixture components for the bulk, the number of iterations of the algorithm, the burn-in period (by default equal to zero) and the thinning interval (by default equal to one).
#' To run the algorithm it is also needed the choice of the starting values, the starting values of the proposal variances, and the parameters of the prior distribution. If not provided, these are automatically set as follows:
#' \itemize{
#' \item \emph{starting values}: \eqn{u} is chosen by the function \code{ithresh} of the \code{threshr} package; \eqn{\xi} and \eqn{\sigma} are chosen by the \code{fpot} function of \code{evd} for data over the threshold; \eqn{\mu} and \eqn{\eta} are chosen as estimates of the \code{gammamixEM} function from the \code{mixtools} package; \eqn{w} is chosen as the vector with entries \eqn{1/k}.
#' \item \emph{variances}: variances for \eqn{\sigma} and \eqn{u} are chosen as the standard deviation of the normal distribution whose 0.01 quantile is equal to 0.9 times the starting value of the associated parameter. The parameters \eqn{\mu_i} and \eqn{\eta_i} are sampled jointly and the proposal variance is chosen using the same method as for \eqn{\sigma} with respect to the parameter \eqn{\mu}. The proposal variance for \eqn{w} is 0.1 and the proposal variance for \eqn{\xi} is 0.1 if negative and 0.25 if positive.
#' \item \emph{prior distributions}: the prior distribution for \eqn{\xi} and \eqn{\sigma} is the objective prior \deqn{p(\xi,\sigma) = \sigma^{-1}(1+\xi)^{-1}(1+2\xi)^{-1/2}.} The prior for the threshold \eqn{u} is Normal with mean chosen as for the starting values and the standard deviation is chosen so that the 0.05 quantile of the prior is equal to the median of the data. The priors for the parameters \eqn{\mu_i} and \eqn{\eta_i} are Gammas with mean chosen as for the starting values and shapes equal to 0.001 so to give high prior variances. The prior for the weigths is the non-informative Dirichlet with parameter 1.
#' }
#'
#' The user can also select any of the three inputs above. \itemize{
#' \item The starting values \code{start} must be a list with entries \code{xi}, \code{sigma}, \code{u}, \code{mu}, \code{eta} and \code{w}. The length of  \code{mu}, \code{eta} and \code{w} must be \code{k}.
#' \item The variances \code{var} must be a list with entries \code{xi}, \code{sigma}, \code{u}, \code{mu} and  \code{w}. The length of \code{mu} must be \code{k}.
#' \item The prior \code{prior} must be a list with entries \code{u}, \code{mu_mu}, \code{mu_eta}, \code{eta_mu} and \code{eta_eta}.  \code{u} gives the mean and the standard deviation of the Normal prior for \eqn{u}. The vectors of length \code{k}, \code{mu_mu} and \code{eta_mu} give the prior means of \eqn{\mu} and \eqn{\eta}, whilst \code{mu_eta} and \code{eta_eta} give their prior shapes.}
#'
#' @references Behrens, Cibele N., Hedibert F. Lopes, and Dani Gamerman. "Bayesian analysis of extreme events with threshold estimation." Statistical Modelling 4.3 (2004): 227-244.
#' @references do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
#'
#' @return \code{fmgpd} returns a list with three elements: \itemize{
#' \item \code{chain}: a matrix of size (\code{it} - \code{burn})/\code{thin}\eqn{\times}5, reporting in each column the posterior sample for each parameter.
#' \item \code{var}: a matrix of size \code{it}\eqn{\times}5 reporting the variance of the proposal distribution for each parameter.
#' \item \code{data}: the dataset used for estimation.
#' }
#' @seealso \code{\link{fggpd}}, \code{\link{mgpd}}
#'
#' @examples \donttest{
#' data(rainfall)
#' ## Small number of iterations and burn-in for quick execution
#' model1 <- fmgpd(rainfall, k = 2, it = 250, burn = 50, thin = 25)
#' start <- list(xi = 0.2, sigma = 2, u = 10, mu = c(2,5), eta = c(2,2) , w = c(0.4,0.6))
#' var <- list(xi = 0.01, sigma = 1, u = 3, mu = c(3,3), w = 0.01)
#' prior <- list(u = c(22,5), mu_mu = c(2,5), mu_eta = c(0.01,0.01),
#'          eta_mu = c(3,3),eta_eta = c(0.01,0.01))
#'
#' model2 <- fmgpd(rainfall, k= 2, it = 250, start = start, var =var, prior = prior)
#' }
#'
#' @param x A vector of positive observations.
#' @param it Number of iterations of the algorithm.
#' @param k number of mixture components for the bulk. Must be either 2, 3, or 4.
#' @param start A list of starting parameter values.
#' @param var A list of starting proposal variance.
#' @param prior A list of hyperparameters for the prior distribution.
#' @param thin Thinning interval.
#' @param burn Burn-in.
#'
#' @export
fmgpd <- function(x, it, k, start = NULL,  var = NULL, prior = NULL, thin = 1, burn = 0){
  if(is.null(start)) start <- starting.values(x,k)
  if(is.null(var)) var <- compute.var(x,start,k)
  if(is.null(prior)) prior <- compute.prior(x,start)

  if(is.null(start) & is.null(prior) & is.null(var)){}
  else{check.input(x, it, k, var, start, prior, thin, burn)}

  start <- unname(unlist(start))
  var <- unname(unlist(var))
  prior <- unname(unlist(prior))
  mu_mu <- prior[3:(3+k-1)]
  mu_eta <- prior[(3+k):(3+2*k-1)]
  eta_mu <- prior[(3+2*k):(3+3*k-1)]
  eta_eta <- prior[(3+3*k):length(prior)]
  prior_u <- prior[1:2]
  prior_mu <-c(rbind(mu_mu, mu_eta))
  prior_eta <- c(rbind(eta_mu, eta_eta))
  start_gpd <- start[1:3]
  start_mu <- start[4:(4+k-1)]
  start_eta <- start[(4+k):(4+2*k-1)]
  start_w <- start[(4+2*k):length(start)]
  mh <- c_fmgpd(x,it,k,start_gpd,start_mu,start_eta,start_w,var,prior_u, prior_mu,prior_eta)
  mh$chain <- mh$chain[seq(burn+1,nrow(mh$chain),by = thin),]
  mh$data <- x
  class(mh) <- c("mgpd","evmm")
  return(mh)
}
