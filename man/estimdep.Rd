\name{estimdep}
\alias{estimdep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dependence estimation
}
\description{
From a set of observations, builds a description of the multivariate distribution
}
\usage{
estimdep(dataframe,varnames,subsampsize,nbsafe=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
a data frame containing the observations
}
  \item{varnames}{
the name of the variables we want to estimate the multivariate distribution
}
  \item{subsampsize}{
the sub-sample size
}
  \item{nbsafe}{
the ratio between the discretized copula size and the number of sub-samples
}
}
\value{
the description of the dependence, it is an object with the following parts:
  \item{cop}{the array representing the discretized copula}
  \item{margins}{the matrix representing the margins, estimated using kernel density estimation}
  \item{varnames}{the names of the variables}
  \item{subsampsize}{as in the function call}
  \item{nbsafe}{as in the function call}
}
\author{
Jerome Collet
}
\examples{
lon=3000
plon=3000
subsampsize=20

##############
x=(runif(lon)-1/2)*3
y=x^2+rnorm(lon)
z=rnorm(lon)
donori=as.data.frame(cbind(x,y,z))
depori=estimdep(donori,c("x","y","z"),subsampsize)

knownvalues=data.frame(z=rnorm(plon))
prev <- predictdep(knownvalues,depori)
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)

knownvalues=data.frame(x=(runif(lon)-1/2)*3)
prev <- predictdep(knownvalues,depori)
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)

knownvalues=data.frame(y=runif(plon,min=-2,max=4))
prev <- predictdep(knownvalues,depori)
plot(prev$x,prev$y,xlim=c(-2,2),ylim=c(-2,5),pch=20,cex=0.5)
points(donori[,1:2],col='red',pch=20,cex=.5)
}
