\name{desscop}
\alias{desscop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete copula graph, a two-dimensional projection
}
\description{
Draws a discrete joint probability, for 2 variables, using bubbles
}
\usage{
desscop(copest, xname, yname, normalize = FALSE, axes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copest}{
the estimated copula (the whole structure resulting from \code{\link{corc}})
}
  \item{xname}{
the name of the variable we want to put on the horizontal axis
}
  \item{yname}{
the name of the variable we want to put on the vertical axis
}
  \item{normalize}{
if TRUE, the smallest probability is rescaled to 0, and the largest to 1
}
  \item{axes}{
if TRUE, puts the name of the variables on the axes
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
tablo = as.data.frame(cbind(x,y))
c=corc(tablo,c("x","y"),8)
desscop(c,"x","y")

tablo = as.data.frame(cbind(x=rep(0,each=lon),y=rep(0,each=lon)))
c=corc(tablo,c("x","y"),8,mixties=TRUE)
desscop(c,"x","y")

}
