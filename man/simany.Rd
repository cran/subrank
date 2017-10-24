\name{simany}
\alias{simany}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test statistic distribution under any hypothesis
}
\description{
Simulates the test statistic, under independence
}
\usage{
simany(sampsize,dimension,subsampsizes,sampnum,nbsafe=5, fun=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampsize}{
sample size
}
  \item{dimension}{
sample dimension
}
  \item{subsampsizes}{
vector of sub-sample sizes
}
  \item{sampnum}{
number of samples
}
  \item{nbsafe}{
the ratio between the number of sub-samples and the cardinality of the discretized copula.
}
  \item{fun}{
the function describing the dependence.
}
  \item{...}{
optional arguments to \code{fun}
}
}
\value{
\item{lrs}{the distances with independent case}
\item{lrs2mean}{the distances with theoretical value, given dependence \code{fun}}
\item{scarcities}{the proportions of non-reached vector ranks}
\item{DistTypes}{a recall of the list of the distance types: "KL","L2","L1","APE"}
}
\author{
Jerome Collet
}
\examples{
depquad <- function(lon,dd,a)
{
  x <- rnorm(lon)
  y0 <- a*x^2
  y <- y0 + rnorm(lon)
  reste=rnorm((dd-2)*lon)
  return(c(x,y,reste))
}
sims0=simany(101,3,8,50,nbsafe=1)
seuils=apply(sims0$lrs,3,quantile,0.95)
seuils=matrix(ncol=4,nrow=50,seuils,byrow=TRUE)
sims1=simany(101,3,8,50,nbsafe=1,fun=depquad,a=0.5)
apply(sims1$lrs[,1,]>seuils,2,mean)
}
