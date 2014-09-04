estimdep <- function(dataframe,varnames,subsampsize,nbsafe=5)
{
  dimension=length(varnames)
  copcomp=corc(dataframe,varnames,subsampsize,nbsafe)
  FdRinv <- list()
  FdR <- list()
  for (var in varnames)
  {
    numvar=pmatch(var,varnames)
    values=sort(unlist(dataframe[var]))
    steps=seq(0,1,along.with=values)
    FdRinv[[numvar]]=approxfun(steps,values,rule=2:2)
    FdR[[numvar]]=approxfun(values,steps,rule=2:2)
  }
  if(copcomp$ties>0.5*copcomp$nsubsampreal)
    warning("To many ties")
  return(list(cop=copcomp$cop,FdRinv=FdRinv,FdR=FdR,
              varnames=varnames,
              subsampsize=subsampsize,nbsafe=nbsafe))
} 

