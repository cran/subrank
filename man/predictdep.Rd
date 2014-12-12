\name{predictdep}
\alias{predictdep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability forecasting
}
\description{
From a set of incomplete observations, and a description of the dependence, provides simulated values of the unknown coordinates. It is also possible to simulate unconditionally, with empty observations.
}
\usage{
predictdep(knownvalues,dependence,smoothing=c("Uniform","Beta"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knownvalues}{
in case of conditional simulation, a matrix containing incomplete observations, the known coordinates being the same for all observations. If no variable name in \code{knwonvalues} appears in \code{dependence$varnames}, then the simulation is unconditional.
}
  \item{dependence}{
the description of the dependence we want to use to forecast, as built by function \code{\link{estimdep}}
}
  \item{smoothing}{
the smoothing method for input and output ranks.
}
}
\value{
the matrix of the completed observations
}
\author{
Jerome Collet
}
\examples{
lon=100
plon=100
subsampsize=10

shift=0
noise=0
knowndims=1

x=rnorm(lon)
y=2*x+noise*rnorm(lon)
donori=as.data.frame(cbind(x,y))
depori=estimdep(donori,c("x","y"),subsampsize)
##
knownvalues=data.frame(x=rnorm(plon)+shift)
prev <- predictdep(knownvalues,depori)
##
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)
##
knownvalues=data.frame(x=rnorm(plon)+shift)
prev <- predictdep(knownvalues,depori,smoothing="Beta")
##
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)

# souci normal si |shift|>>1

knownvalues=data.frame(z=rnorm(plon)+shift)
prev <- predictdep(knownvalues,depori)
##
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)
}
