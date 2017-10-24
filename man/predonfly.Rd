\name{predonfly}
\alias{predonfly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability forecasting
}
\description{
From two sets of observations, first one of complete observations and second one of incomplete observations, provides simulated values of the unknown coordinates.
}
\usage{
predonfly(completeobs,incompleteobs,varnames,subsampsize,nbpreds=1,mixties=FALSE,
          maxtirs=1e5,complete=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{completeobs}{
the set of complete observations.
}
  \item{incompleteobs}{
the set of incomplete observations.
}
  \item{varnames}{
the modeled variables.
}
  \item{subsampsize}{
the sub-sample size.
}
  \item{nbpreds}{
the number of predictions for each incomplete observation.
}
  \item{mixties}{
if \code{TRUE}, should put equal weight on tied values, using random permutations (not yet implemented)
}
  \item{maxtirs}{
the maximum number of sub-samples, to stop the computation even if they did not provide \code{nbpreds} predictions for each incomplete observation.
}
  \item{complete}{
If \code{TRUE}, predictions are completed with incomplete observations
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
plon=30
subsampsize=10

x=rnorm(lon)
y=2*x+rnorm(lon)*0
donori=as.data.frame(cbind(x,y))
##
knownvalues=data.frame(x=rnorm(plon))
prev <- predonfly(donori,knownvalues,c("x","y"),subsampsize,100)

##
plot(prev$x,prev$y,pch=20,cex=0.5,
     ylim=range(c(prev$y,donori$y),na.rm=TRUE),xlim=range(c(prev$x,donori$x)))
points(donori[,1:2],col='red',pch=20,cex=.5)

lon=3000
mg=20
dimtot=4
rayon=6

genboules <- function(lon,a,d)
{
  ss <- function(vec)
  {return(sum(vec*vec))}
  surface=matrix(nrow=lon,ncol=d,data=rnorm(lon*d))
  rayons=sqrt(apply(surface,1,ss))
  surface=surface/rayons
  return(matrix(nrow=lon,ncol=d,data=rnorm(lon*d))+a*surface)
}

##############

donori=genboules(lon,rayon,dimtot)
donori=as.data.frame(donori)

dimconnues=3:dimtot
valconnues=matrix(nrow=1,ncol=length(dimconnues),data=0)
valconnues=as.data.frame(valconnues)
names(valconnues)=names(donori)[3:dimtot]
prev <- predonfly(donori,valconnues,names(donori),subsampsize,100)

boule2=genboules(plon,rayon,2)

plot(boule2[,1:2],xlab='X1',ylab='X2',pch=20,cex=.5)
plot(prev$V1,prev$V2,xlab='X1',ylab='X2',pch=20,cex=.5)
}
