\name{desscoptous}
\alias{desscoptous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete copula graph, ALL two-dimensional projections
}
\description{
Draws a discrete joint probability, for 2 variables, using bubbles
}
\usage{
desscoptous(copest, normalize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copest}{
the estimated copula (the whole structure resulting from \code{\link{corc}})
}
  \item{normalize}{
if TRUE, the smallest probability is rescaled to 0, and the largest to 1
}
}
\author{
Jerome Collet
}
\examples{
lon <- 31
a <- 2.85
x <- rnorm(lon)
y = a*x^2+rnorm(lon)
z = rnorm(lon)
tablo = as.data.frame(cbind(x,y,z))
c=corc(tablo,c("x","y","z"),8)
desscoptous(c)
}
