\name{subrank-package}
\alias{subrank-package}
\alias{subrank}
\docType{package}
\title{
Package to estimate copula using ranks and sub-sampling
}
\description{
Taking a sample, its dimension, and a sub-sample size, allows to estimate a discretized copula. This object has interesting features: convergence to copula, robustness with respect to dimension.
}
\details{
\tabular{ll}{
Package: \tab subrank\cr
Type: \tab Package\cr
Version: \tab 0.9.4\cr
Date: \tab 2014-11-11\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}
The most important functions are corc and predonfly: first one estimates the discrete copula, second predicts using a copula estimated on the fly..
}
\author{
Jerome Collet

Maintainer: Jerome Collet <jeromepcollet@gmail.com>
}
\keyword{ package }
\examples{
lon <- 31
a <- 2.85
x <- rnorm(lon)
y = a*x^2+rnorm(lon)
tablo = as.data.frame(cbind(x,y))
c=corc(tablo,c(1,2),8)
desscop(c,1,2)
}
