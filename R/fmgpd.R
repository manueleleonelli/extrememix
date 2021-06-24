#' @import mixtools 
#' @import evd
#' @importFrom stats quantile
#' @importFrom utils capture.output
starting.values <- function(x, k){
  lower <- subset(x, x <= quantile(x,0.9))
  pargpd <- unname(fpot(x,quantile(x,0.9),std.err = F)$estimate)
  invisible(capture.output(mix <- gammamixEM(x, k = k, epsilon = 1e-02)))
  mu <- unname(mix$gamma.pars[1,])*unname(mix$gamma.pars[2,])
  eta <- unname(mix$gamma.pars[1,])
  return(list(startxi = max(-0.4,pargpd[2]),startsigma = pargpd[1], startu = unname(quantile(x,0.9)), startmu = sort(mu) , starteta = eta[order(mu)], startw = mix$lambda[order(mu)]))
}


#' @importFrom stats uniroot
#' @importFrom stats qnorm
compute.var <- function(x,start,k){
  if(start$startxi>=0 ){varxi <- 0.01} else{varxi <- 0.001}
  varsigma <- uniroot(function(r) qnorm(0.01, start$startsigma,r) - max(0.9*start$startsigma,0) , c(0,max(x)))$root
  varu <- uniroot(function(r) qnorm(0.01,start$startu,r) - max(0.9*start$startu,0), c(0,max(x)))$root
  varmu <- rep(0,length(start$startmu))
  for(i in 1:k){
  varmu[i] <- uniroot(function(r) qnorm(0.01,start$startmu[i],r) -  max(0.99*start$startmu[i],0), c(0,max(x)))$root
  }
  varw <- 10
  return(list(Vxi=varxi, Vsigma = varsigma, Vu = varu, Vmu = varmu, Vw = varw))
}


check.input <- function(x, it, k, var, start, prior, thin, burn){
  if(any(x <=0)) stop("x must be positive")
  if(k%% 1 !=0 | k <= 0 ) stop("k must be a positive integer")
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
  if(is.null(var[["Vw"]])){stop('Vw not provided')}
 
  if(is.null(start[["startxi"]])){stop('startxi not provided')}
  if(is.null(start[["startsigma"]])){stop('startsigma not provided')}
  if(is.null(start[["startu"]])){stop('startu not provided')}
  if(is.null(start[["startmu"]])){stop('startmu not provided')}
  if(is.null(start[["starteta"]])){stop('starteta not provided')}
  if(is.null(start[["startw"]])){stop('startw not provided')}

  if(is.null(prior[["prior_u"]])){stop('Prior for u not provided')}
  if(is.null(prior[["mu_mu"]])){stop('Prior mean for mu not provided')}
  if(is.null(prior[["mu_eta"]])){stop('Prior shape for mu not provided')}
  if(is.null(prior[["eta_mu"]])){stop('Prior mean for eta not provided')}
  if(is.null(prior[["eta_eta"]])){stop('Prior shape for eta not provided')}
}


#' @importFrom stats uniroot
#' @importFrom stats qnorm
compute.prior <- function(x,start){
  eta_mu <- start$starteta 
  eta_eta <- rep(0.001, length(start$starteta))
  prior_u <- c(start$startu,  uniroot(function(r) qnorm(0.025,start$startu,r) - quantile(x,0.5), c(0,max(x)))$root)
  mu_mu <- start$startmu
  mu_eta <- rep(0.001, length(start$starteta))
  return(list(prior_u = prior_u, mu_mu = mu_mu, mu_eta = mu_eta, eta_mu = eta_mu, eta_eta = eta_eta))
}


#' Fit of the GGPD model using an MCMC algorithm.
#' @param x A vector of positive observations.
#' @param it Number of iterations of the algorithm.
#' @param k number of mixture components for the bulk
#' @param start A list of starting parameter values.
#' @param var A list of starting proposal variance.
#' @param prior A list of hyperparameters for the prior distribution.
#' @param thin Thinning interval.
#' @param burn Burn-in.
#' 
#' @export
#' @export
fmgpd <- function(x, it, k, start = NULL,  var = NULL, prior = NULL, thin = 1, burn = 0){
  if(k == 1) warning("use of fggpd is recommended for k = 1")
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
  class(mh) <- "mgpd"
 return(mh)
}
