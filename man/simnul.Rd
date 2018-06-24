\name{simnul}
\alias{simnul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test statistic distribution under independence hypothesis
}
\description{
Simulates the test statistic, under independence
}
\usage{
simnul(sampsize, dimension, subsampsizes, sampnum,KL=TRUE,nbsafe=5,nthreads=2)
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
  \item{KL}{
if TRUE, returns the Kullback-Leibler divergence with the independent case, if FALSE, the L2 distance. There is no re-normalization, contrary to what happens for \code{simany}.
}
  \item{nbsafe}{
the ratio between the number of sub-samples and the cardinality of the discretized copula.
}
  \item{nthreads}{
number of number of threads, assumed to be strictly positive. For "full throttle" computations, consider using parallel::detectCores()}
}
\value{
\item{lrs }{the distances with independent case}
\item{scarcities }{the proportions of non-reached vector ranks}
}
\author{
Jerome Collet
}
\examples{
library(datasets)
# plot(swiss)
c=corc(swiss,1:3,8)
c
RV=sum(c$cop*log(c$cop),na.rm=TRUE)+3*log(8)
sims=simnul(47,3,8,100)
pvalue=mean(RV<sims$lrs)
pvalue
RV
summary(sims$lrs)
}
