corc0 <-
function(datavector,sampsize,dimension,subsampsize,nboot,u,mixties=FALSE,nthreads=2)
{
  return( .Call("InterCopulation", PACKAGE = "subrank",
    as.double(datavector),
    as.integer(sampsize), as.integer(subsampsize), as.integer(dimension),
    as.integer(mixties),
    as.integer(u), as.integer(nboot),
	as.integer(nthreads)
               )
        )
}
