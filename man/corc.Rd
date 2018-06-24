\name{corc}
\alias{corc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate copula using ranks and sub-sampling
}
\description{
Takes a sample, its dimension, a sub-sample size, and returns a discrete copula.
}
\usage{
corc(dataframe, varnames, subsampsize, nbsafe=5,mixties=FALSE,nthreads=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
a data frame, containing the observations
}
  \item{varnames}{
the name of the variables we want to estimate the dependence between
}
  \item{subsampsize}{
the sub-sample size
}
  \item{nbsafe}{
the ratio between the number of sub-samples and the cardinality of the discretized copula.
}
  \item{mixties}{
if \code{TRUE}, put equal weight on tied values, using random permutations 
}
  \item{nthreads}{
number of number of threads, assumed to be strictly positive. For "full throttle" computations, consider using parallel::detectCores()}
}
\value{
\item{cop}{an array representing the discretized copula}
\item{ties}{the number of sub-samples with a tie}
\item{nsubsampreal}{the effective number of sub-samples drawn}
\item{varnames}{the name of the variables studied}
\item{nnm}{the number of observations without missing values}
}
\author{
Jerome Collet
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
lon <- 30
a <- 2
x <- rnorm(lon)
y = a*x^2+rnorm(lon)
datatable = as.data.frame(cbind(x,y))
c=corc(datatable,c("x","y"),8)
c
sum(c$cop)
}
