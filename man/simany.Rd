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
simany(sampsize,subsampsizes,sampnum,nbsafe=5, fun=5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampsize}{
sample size
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
if numeric, indicates independence in dimension \code{fun}, elsewhere the function describing the dependence.
}
  \item{...}{
optional arguments to \code{fun}
}
}
\value{
\item{lrs}{the distances with independent case}
\item{lrs2mean}{the distances with asymptotic value, given dependence \code{fun}}
\item{lrsmeanmean}{the distances between asymptotic values, given dependence \code{fun} and in independent case}
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
  return(matrix(ncol=dd,c(x,y,reste)))
}
sims0=simany(101,5,50,nbsafe=1,fun=3)
seuils=apply(sims0$lrs,3,quantile,0.9)
seuils=matrix(ncol=4,nrow=10,seuils,byrow=TRUE)
sims1=simany(101,5,10,nbsafe=1,fun=depquad,a=0.3,dd=3)
apply(sims1$lrs[,1,]>seuils,2,mean)
}
