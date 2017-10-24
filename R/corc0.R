corc0 <-
function(datavector,sampsize,dimension,subsampsize,nboot,u,mixties=FALSE)
{
  return( .Call("InterCopulation",
    as.double(datavector),
    as.integer(sampsize), as.integer(subsampsize), as.integer(dimension),
    as.integer(mixties),
    as.integer(u), as.integer(nboot) 
               )
        )

}
