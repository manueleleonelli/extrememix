% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggpd.R
\name{ggpd}
\alias{ggpd}
\alias{dggpd}
\alias{pggpd}
\alias{qggpd}
\alias{rggpd}
\title{The GGPD distribution}
\usage{
dggpd(x, xi, sigma, u, mu, eta, log = FALSE)

pggpd(q, xi, sigma, u, mu, eta, lower.tail = TRUE)

qggpd(p, xi, sigma, u, mu, eta, lower.tail = TRUE)

rggpd(N, xi, sigma, u, mu, eta)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{xi}{shape parameter of the tail GPD (scalar).}

\item{sigma}{scale parameter of the tail GPD (scalar).}

\item{u}{threshold parameter of the tail GPD (scalar).}

\item{mu}{mean of the gamma bulk (scalar).}

\item{eta}{shape of the gamma bulk (scalar).}

\item{log}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.}

\item{p}{vector of probabilities.}

\item{N}{number of observations.}
}
\value{
The GGPD distribution is an extreme value mixture model with density
\deqn{f_{GGPD}(x|\xi,\sigma,u,\mu,\eta,w)=\left\{\begin{array}{ll} f_{GA}(x|\mu,\eta), & x\leq u \\ (1-F_{GA}(u|\mu,\eta))f_{GPD}(x|\xi,\sigma,u), &\mbox{otherwise},  \end{array}\right.} where \eqn{f_{GA}} is the density of the Gamma parametrized by mean \eqn{\mu} and shape \eqn{\eta}, \eqn{F_{GA}} is the distribution function of the Gamma and \eqn{f_{GPD}} is the density of the Generalized Pareto Distribution, i.e.
\deqn{f_{GPD}(x|\xi,\sigma,u)=\left\{\begin{array}{ll} 1- (1+\frac{\xi}{\sigma}(x-u))^{-1/\xi}, & \mbox{if } \xi\neq 0,\\ 1- \exp\left(-\frac{x-u}{\sigma}\right), & \mbox{if } \xi = 0, \end{array}\right.}
where \eqn{\xi} is a shape parameter, \eqn{\sigma > 0} is a scale parameter and \eqn{u>0} is a threshold.

\code{dggpd} gives the density, \code{pggpd} gives the distribution function, \code{qggpd} gives the quantile function, and \code{rggpd} generates random deviates. The length of the result is determined by \code{N} for \code{rggpd} and by the length of \code{x}, \code{q} or \code{p} otherwise.
}
\description{
Density, distribution function, quantile function and random generation for the GGPD distribution.
}
\examples{
dggpd(3, xi = 0.5, sigma = 2, u = 5, mu = 3, eta = 3)


}
\references{
Behrens, Cibele N., Hedibert F. Lopes, and Dani Gamerman. "Bayesian analysis of extreme events with threshold estimation." Statistical Modelling 4.3 (2004): 227-244.
}
