\name{subrank-package}
\alias{subrank-package}
\alias{subrank}
\docType{package}
\title{
\packageTitle{subrank}
}
\description{
\packageDescription{subrank}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{subrank}
\packageIndices{subrank}
Taking a sample, its dimension, and a sub-sample size, allows to estimate a discretized copula. This object has interesting features: convergence to copula, robustness with respect to dimension.
}
\author{
\packageAuthor{subrank}

Maintainer: \packageMaintainer{subrank}
}
\examples{
lon <- 31
a <- 2.85
x <- rnorm(lon)
y = a*x^2+rnorm(lon)
tablo = as.data.frame(cbind(x,y))
c=corc(tablo,c(1,2),8)
desscop(c,1,2)
}
