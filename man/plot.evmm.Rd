% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.evmm}
\alias{plot.evmm}
\title{Plot of Extreme Value Mixture Models}
\usage{
\method{plot}{evmm}(x, ...)
}
\arguments{
\item{x}{an object of class \code{evmm}.}

\item{...}{additional parameters for compatibility.}
}
\description{
Plotting method for objects of class \code{evmm} giving an overview of an estimated model.
}
\details{
The \code{plot} method for objects of class \code{evmm} reports four plots: \itemize{
\item An histogram of the posterior distribution of xi.
\item An histogram of the posterior distribution of sigma.
\item A line plot of the estimated quantiles, their posterior credibility interval, and the empirical ones.
\item A plot of the predictive distribution together with the data histogram.
}
}
\examples{
\dontrun{
data(rainfall)
model1 <- fggpd(rainfall, it = 25000, burn = 5000, thin = 25)
plot(model1)
}

}
