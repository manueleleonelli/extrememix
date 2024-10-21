#' @importFrom threshr ithresh
#' @importFrom stats var
#' @importFrom evd fpot
ggpd.starting.values <- function(x){
  u <- unname(summary(ithresh(data = x, u_vec =  stats::quantile(x, probs = seq(0.75, 0.95, by = 0.05)), n = 100)))[3]
  lower <- subset(x, x <= u)
  pargpd <- unname(fpot(x,u,std.err = F)$estimate)
  return(list(xi = max(-0.4,pargpd[2]),sigma = pargpd[1], u = u, mu = mean(lower), eta = mean(lower)^2/var(lower)))
}

#' @importFrom stats qnorm
#' @importFrom stats uniroot
ggpd.compute.var <- function(x,start){
  if(start$xi>=0 ){xi <- 0.025} else{xi <- 0.001}
  sigma <- uniroot(function(r) qnorm(0.01, start$sigma,r) - max(0.9*start$sigma,0) , c(0,max(x)))$root
  u <- uniroot(function(r) qnorm(0.01,start$u,r) - max(0.9*start$u,0), c(0,max(x)))$root
  eta <- uniroot(function(r) qnorm(0.01,start$mu,r)- max(0.99*start$mu,0), c(0,max(x)))$root
  mu <- uniroot(function(r) qnorm(0.01,start$eta,r) -  max(0.99*start$eta,0), c(0,max(x)))$root
  return(list(xi=xi, sigma = sigma, u = u, mu = mu, eta = eta))
}

#' @importFrom stats quantile
#' @importFrom stats uniroot
ggpd.compute.prior <- function(x,start){
  eta <- c(start$eta, 0.001)
  u <- c(start$u, uniroot(function(r) qnorm(0.05,start$u,r) - stats::quantile(x,0.5), c(0,max(x)))$root)
  mu <- c(start$mu, 0.001)
  return(list(u = u, mu = mu, eta = eta))
}


ggpd.check.input <- function(it,var,start,prior, thin, burn){
  if(it %% 1 != 0 | it <= 0){stop('iterations must be a positive integer')}
  if(burn > it){stop('burn-in cannot be larger than iterations')}
  if(thin > (it-burn)){stop('thinning too large')}
  if(mode(var) != "list"){stop('var must be a list')}
  if(mode(start) != "list"){stop('start must be a list')}
  if(mode(prior) != "list"){stop('prior must be a list')}

  if(is.null(var[["xi"]])){stop('Variance for xi not provided')}
  if(is.null(var[["sigma"]])){stop('Variance for sigma not provided')}
  if(is.null(var[["u"]])){stop('Variance for u not provided')}
  if(is.null(var[["mu"]])){stop('Variance for mu not provided')}
  if(is.null(var[["eta"]])){stop('Variance for eta not provided')}
  if(any(var<0)){stop('Proposal variances must be positive')}

  if(is.null(start[["xi"]])){stop('starting value for xi not provided')}
  if(is.null(start[["sigma"]])){stop('starting value for sigma not provided')}
  if(is.null(start[["u"]])){stop('starting value for u not provided')}
  if(is.null(start[["mu"]])){stop('starting value for mu not provided')}
  if(is.null(start[["eta"]])){stop('starting value for eta not provided')}
  if(any(start[-1]<0)){stop('Starting values for sigma, u, mu and eta must be positive')}
  if(start[1]< -0.5){stop('Starting value for xi must be larger than -0.5')}

  if(is.null(prior[["u"]])){stop('Prior for u not provided')}
  if(is.null(prior[["mu"]])){stop('Prior for mu not provided')}
  if(is.null(prior[["eta"]])){stop('Prior for eta not provided')}
  if(length(prior[["u"]])!= 2) stop("Two parameters must be given to u prior")
  if(length(prior[["mu"]])!= 2) stop("Two parameters must be given to mu prior")
  if(length(prior[["eta"]])!= 2) stop("Two parameters must be given to eta prior")
  if(any(prior[["u"]] <0) || any(prior[["mu"]] <0) || any(prior[["eta"]]<0) ){stop('Prior hyperparameters must be positive')}
}

#' GGPD Estimation
#'
#' Fit of the GGPD model using an MCMC algorithm.
#'
#' Estimation of the GGPD is carried out using an adaptive block Metropolis-Hastings algorithm. As standard, the user needs to specify the data to use during estimation, the number of iterations of the algorithm, the burn-in period (by default equal to zero) and the thinning interval (by default equal to one).
#' To run the algorithm it is also needed the choice of the starting values, the starting values of the proposal variances, and the parameters of the prior distribution. If not provided, these are automatically set as follows:
#' \itemize{
#' \item \emph{starting values}: \eqn{u} is chosen by the function \code{ithresh} of the \code{threshr} package; \eqn{\xi} and \eqn{\sigma} are chosen by the \code{fpot} function of \code{evd} for data over the threshold; \eqn{\mu} and \eqn{\eta} are chosen as the maximum likelihood estimate of the Gamma distribution over data below the threshold.
#' \item \emph{variances}: variances are chosen as the standard deviation of the normal distribution whose 0.01 quantile is equal to 0.9 times the starting value of the associated parameter.
#' \item \emph{prior distributions}: the prior distribution for \eqn{\xi} and \eqn{\sigma} is the objective prior \deqn{p(\xi,\sigma) = \sigma^{-1}(1+\xi)^{-1}(1+2\xi)^{-1/2}.} The prior for the threshold \eqn{u} is Normal with mean chosen as for the starting values and the standard deviation is chosen so that the 0.05 quantile of the prior is equal to the median of the data. The priors for the parameters \eqn{\mu} and \eqn{\eta} are Gammas with mean chosen as for the starting values and shapes equal to 0.001 so to give high prior variances.
#' }
#'
#' The user can also select any of the three inputs above. \itemize{
#' \item The starting values \code{start} must be a list with entries \code{xi}, \code{sigma}, \code{u}, \code{mu}, \code{eta}.
#' \item The variances \code{var} must be a list with entries \code{xi}, \code{sigma}, \code{u}, \code{mu}, \code{eta}.
#' \item The prior \code{prior} must be a list with entries \code{u}, \code{mu}, \code{eta} all containing a vector of length two (for \code{u} giving the mean and the standard deviation of the Normal prior, for \code{mu} and \code{eta} giving the mean and shape of the Gamma prior).}
#'
#' @examples \donttest{
#' ## Small number of iterations and burn-in for quick execution
#' data(rainfall)
#' model1 <- fggpd(rainfall, it = 250, burn = 50, thin = 25)
#'
#' start <- list(xi = 0.2, sigma = 2, u = 10, mu = 5, eta = 2)
#' var <- list(xi = 0.01, sigma = 1, u = 3, mu = 3, eta = 1)
#' prior <- list(u = c(22,5), mu = c(4,16), eta = c(0.001,0.001))
#' model2 <- fggpd(rainfall,it = 250, start = start, var =var, prior = prior)
#' }
#'
#' @references Behrens, Cibele N., Hedibert F. Lopes, and Dani Gamerman. "Bayesian analysis of extreme events with threshold estimation." Statistical Modelling 4.3 (2004): 227-244.
#' @references do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
#'
#' @return \code{fggpd} returns a list with three elements: \itemize{
#' \item \code{chain}: a matrix of size (\code{it} - \code{burn})/\code{thin}\eqn{\times}5, reporting in each column the posterior sample for each parameter.
#' \item \code{var}: a matrix of size \code{it}\eqn{\times}5 reporting the variance of the proposal distribution for each parameter.
#' \item \code{data}: the dataset used for estimation.
#' }
#'
#' @seealso \code{\link{ggpd}}
#' @param x A vector of positive observations.
#' @param it Number of iterations of the algorithm.
#' @param start A list of starting parameter values.
#' @param var A list of starting proposal variances.
#' @param prior A list of hyperparameters for the prior distribution.
#' @param thin Thinning interval.
#' @param burn Burn-in length.
#' @import RcppProgress
#' @export
fggpd <- function(x, it,start = NULL, var = NULL,  prior = NULL, thin = 1, burn = 0){
  if(is.null(start)) start <- ggpd.starting.values(x)
  if(is.null(var)) var <- ggpd.compute.var(x,start)
  if(is.null(prior)) prior <- ggpd.compute.prior(x,start)
  ggpd.check.input(it,var,start,prior,thin, burn)

  start <- c(start$xi,start$sigma, start$u, start$mu, start$eta)
  var <- c(var$xi, var$sigma, var$u, var$mu, var$eta)
  prior <- c(prior$u, prior$mu, prior$eta)

  mh <- c_fggpd(x,it,start,var,prior)
  mh$chain <- mh$chain[seq(burn+2,nrow(mh$chain),by = thin),]
  colnames(mh$chain) <- c("xi","sigma","u","mu","eta")
  colnames(mh$chain) <- c("xi","sigma","u","mu","eta")
  mh$data <- x
  class(mh) <- c("ggpd","evmm")
  return(mh)
}
