\name{corc0}
\alias{corc0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate copula using ranks and sub-sampling, minimal version.
}
\description{
Minimal version of function corc.
}
\usage{
corc0(datavector,sampsize,dimension,subsampsize,nboot,u,mixties=FALSE,nthreads=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datavector}{
a vector, containing the observations
}
  \item{sampsize}{
the sample size
}
  \item{dimension}{
the sample dimension
}
  \item{subsampsize}{
the sub-sample size
}
  \item{nboot}{
the number of sub-samples (must be big)
}
  \item{u}{
a random seed, integer
}
  \item{mixties}{
if \code{TRUE}, put equal weight on tied values, using random permutations 
}
  \item{nthreads}{
number of number of threads, assumed to be strictly positive. For "full throttle" computations, consider using parallel::detectCores()}
}
\value{
the number of hits for each vector of ranks, plus 2 last values of the vector : number of ties and number of sub-samples really used.
}
\author{
Jerome Collet
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
lon <- 30
a <- 2.85
x <- rnorm(lon)
y = a*x^2+rnorm(lon)
c=corc0(c(x,y),lon,2,8,1e5,75014)
c

c0=c(
1203, 1671, 1766, 959, 1586, 1715, 1803, 1205, 1260,1988, 2348, 1917, 3506, 2045, 1340,
1093, 2694, 2757,2233, 1085, 2322, 1793, 1569, 1263, 1709, 1747, 1512,1308, 1778, 1354,
1184, 1097, 2487, 2730, 2112, 1100,2435, 2033, 1572, 1093, 1369, 1722, 1462, 1015, 1228,
1419, 1776, 1852, 1009, 1097, 1179, 1323, 1595, 1316,1477, 2628, 889, 1178, 1981, 4000, 
35, 840, 2091, 4467,0, 27405)
set.seed(75013)
lon=30
dimension=3
sssize=4
c0==corc0(rnorm(lon*dimension),lon,dimension,sssize,1e5,75014)
}
