\name{extrememix-package}
\alias{extrememix-package}
\alias{extrememix}
\docType{package}
\title{
  Bayesian Estimation of Extreme Value Mixture Models
}
\description{
  Bayesian estimation of various extreme value mixture models using MCMC algorithms. The package provides capabilities for model selection and comparison, estimation of quantities of interest in extreme value analyses and plotting functions taking advantage of the \code{ggplot2} package.
}
\details{
  The functions in \code{extrememix} can be divided into four groups according to their purpose:
  
  \itemize{
  \item Density functions, probability functions, quantile functions and random number generation of many extreme value mixture models:
  \itemize{
  \item Gamma bulk and GPD tail (\code{\link{ggpd}})
  \item Mixture of Gammas bulk and GPD tail (\code{\link{mgpd}})
  }
  \item Bayesian estimation using MCMC algorithms:
  \itemize{
  \item \code{\link{fggpd}}: estimation of the \code{\link{ggpd}} model
  \item \code{\link{fmgpd}}: estimation of the \code{\link{mgpd}} model
  }
  \item Model selection functions:
  \itemize{
  \item \code{AIC}: Aikiake information criterion
  \item \code{BIC}: Bayesian information criterion
  \item \code{\link{DIC}}: Deviance information criterion
  \item \code{\link{WAIC}}: Widely applicable information criterion
  }
  \item Summaries for extreme value analysis:
  \itemize{
  \item \code{\link{pred}}: plot of the predictive distribution
  \item \code{\link{quant}}: estimated and empirical quantiles
  \item \code{\link{return_level}}: estimated and empirical return levels
  \item \code{\link{VaR}}: estimated and empirical value-at-risk
  \item \code{\link{ES}}: estimated and empirical expected shortfall
  }
  }
All estimated models as well as summaries have their own printing and plotting method.
}
\author{
Manuele Leonelli <manuele.leonelli@ie.edu>
}
\references{
  Behrens, Cibele N., Hedibert F. Lopes, and Dani Gamerman. "Bayesian analysis of extreme events with threshold estimation." Statistical Modelling 4.3 (2004): 227-244.
  
  do Nascimento, Fernando Ferraz, Dani Gamerman, and Hedibert Freitas Lopes. "A semiparametric Bayesian approach to extreme value estimation." Statistics and Computing 22.2 (2012): 661-675.
  
    Li, Yunxian, Niansheng Tang, and Xuejun Jiang. "Bayesian approaches for analyzing earthquake catastrophic risk." Insurance: Mathematics and Economics 68 (2016): 110-119.
    
  Scarrott, Carl, and Anna MacDonald. "A review of extreme value threshold es-timation and uncertainty quantification." REVSTAT–Statistical Journal 10.1 (2012): 33-60.

}
\examples{
  \dontrun{
        data(rainfall)
        model <- fggpd(rainfall, it = 25000, burn = 5000, thin = 25)
        plot(model)
        WAIC(model)
        quant(model)
      }
}





