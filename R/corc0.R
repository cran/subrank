corc0 <-
function(datavector,sampsize,dimension,subsampsize,nboot,u)
{
  return( .Call("InterCopulation",
    as.double(datavector),
    as.integer(sampsize), as.integer(subsampsize), as.integer(dimension),
    as.integer(u), as.integer(nboot) 
               )
        )

}
