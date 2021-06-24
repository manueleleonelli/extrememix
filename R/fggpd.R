#' @importFrom stats var
ggpd.starting.values <- function(x){
  lower <- subset(x, x <= quantile(x,0.9))
  pargpd <- unname(fpot(x,quantile(x,0.9),std.err = F)$estimate)
  return(list(startxi = max(-0.4,pargpd[2]),startsigma = pargpd[1], startu = unname(quantile(x,0.9)), startmu = mean(lower), starteta = mean(lower)^2/var(lower)))
}


ggpd.compute.var <- function(x,start){
  varxi <- 0.01
  varsigma <- uniroot(function(r) qnorm(0.01, start$startsigma,r) - max(0.9*start$startsigma,0) , c(0,max(x)))$root
  varu <- uniroot(function(r) qnorm(0.01,start$startu,r) - max(0.9*start$startu,0), c(0,max(x)))$root
  vareta <- uniroot(function(r) qnorm(0.01,start$startmu,r)- max(0.99*start$startmu,0), c(0,max(x)))$root
  varmu <- uniroot(function(r) qnorm(0.01,start$starteta,r) -  max(0.99*start$starteta,0), c(0,max(x)))$root
  return(list(Vxi=varxi, Vsigma = varsigma, Vu = varu, Vmu = varmu, Veta = vareta))
}


ggpd.compute.prior <- function(x,start){
  prior.eta <- list(c = start$starteta, d = 0.001)
  prior.u <- list(mean = start$startu, sd = uniroot(function(r) qnorm(0.025,start$startu,r) - quantile(x,0.5), c(0,max(x)))$root)
  prior.mu <- list(a = start$startmu, b = 0.001)
  return(list(u = prior.u, mu = prior.mu, eta = prior.eta))
}


ggpd.check.input <- function(it,var,start,prior, thin, burn){
  if(it %% 1 != 0 | it <= 0){stop('iterations must be a positive integer')}
  if(burn > it){stop('burn-in cannot be larger than iterations')}
  if(thin > (it-burn)){stop('thinning too large')}
  if(mode(var) != "list"){stop('var must be a list')}
  if(mode(start) != "list"){stop('start must be a list')}
  if(mode(prior) != "list"){stop('prior must be a list')}
  
  if(is.null(var[["Vxi"]])){stop('Vxi not provided')}
  if(is.null(var[["Vsigma"]])){stop('Vsigma not provided')}
  if(is.null(var[["Vu"]])){stop('Vu not provided')}
  if(is.null(var[["Vmu"]])){stop('Vmu not provided')}
  if(is.null(var[["Veta"]])){stop('Veta not provided')}
  if(any(var<0)){stop('Proposal variances must be positive')}
  
  if(is.null(start[["startxi"]])){stop('startxi not provided')}
  if(is.null(start[["startsigma"]])){stop('startsigma not provided')}
  if(is.null(start[["startu"]])){stop('startu not provided')}
  if(is.null(start[["startmu"]])){stop('startmu not provided')}
  if(is.null(start[["starteta"]])){stop('starteta not provided')}
  if(any(start[-1]<0)){stop('Starting values for sigma, u, mu and eta must be positive')}
  if(start[1]< -0.5){stop('Starting value for xi must be larger than -0.5')}
  
  if(is.null(prior$u[["mean"]])){stop('Prior mean for u not provided')}
  if(is.null(prior$u[["sd"]])){stop('Prior sd for u not provided')}
  if(is.null(prior$mu[["a"]])){stop('Prior parameter a for mu not provided')}
  if(is.null(prior$mu[["b"]])){stop('Prior parameter b for mu not provided')}
  if(is.null(prior$eta[["c"]])){stop('Prior parameter c for eta not provided')}
  if(is.null(prior$eta[["d"]])){stop('Prior parameter d for eta not provided')}
  if(any(prior$u <0) || any(prior$mu <0) || any(prior$eta<0) ){stop('Prior hyperparameters must be positive')}
}

#' Fit of the GGPD model using an MCMC algorithm.
#' @param x A vector of positive observations.
#' @param it Number of iterations of the algorithm.
#' @param start A list of starting parameter values.
#' @param var A list of starting proposal variance.
#' @param prior A list of hyperparameters for the prior distribution.
#' @param thin Thinning interval.
#' @param burn Burn-in.
#' 
#' @export
fggpd <- function(x, it,start = NULL, var = NULL,  prior = NULL, thin = 1, burn = 0){
  if(is.null(start)) start <- ggpd.starting.values(x)
  if(is.null(var)) var <- ggpd.compute.var(x,start)
  if(is.null(prior)) prior <- ggpd.compute.prior(x,start)
  ggpd.check.input(it,var,start,prior,thin, burn)
  
  start <- c(start$startxi,start$startsigma, start$startu, start$startmu, start$starteta)
  var <- c(var$Vxi, var$Vsigma, var$Vu, var$Vmu, var$Veta)
  prior <- c(prior$u$mean, prior$u$sd, prior$mu$a, prior$mu$b, prior$eta$c,prior$eta$d)
  
  mh <- c_fggpd(x,it,start,var,prior)
  mh$chain <- mh$chain[seq(burn+1,nrow(mh$chain),by = thin),]
  mh$data <- x
  class(mh) <- "ggpd"
  return(mh)
}