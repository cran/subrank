simany <- function(sampsize,dimension,subsampsizes,sampnum,nbsafe=5, fun=NULL, ...)
{
  DistTypes=colnames(distance(1:5,5:1))
  lrs=array(dim=c(sampnum,length(subsampsizes),length(DistTypes)))
  lrs2mean=array(dim=c(sampnum,length(subsampsizes),length(DistTypes)))
  scarcities=matrix(ncol=length(subsampsizes),nrow=sampnum)
  if (!is.null(fun)) fun <- match.fun(fun)
  vraicop=list()
  for (s in 1:length(subsampsizes))
  {
    nboot=nbsafe*subsampsizes[s]^dimension
    if (is.null(fun)) {
      vraicop[[s]]=rep(1/subsampsizes[s]^dimension,subsampsizes[s]^dimension)
    }
    else {
      simdata=simdata=fun(sampsize*sampnum,dimension, ...)
      vraicoptemp=corc0(simdata,sampsize*sampnum,dimension,subsampsizes[s],nboot,42)
      nbootreel=vraicoptemp[subsampsizes[s]^dimension+2]
      vraicop[[s]]=vraicoptemp[1:subsampsizes[s]^dimension]/nbootreel
    }
  }
  for (e in 1:sampnum)
  {
    for (s in 1:length(subsampsizes))
    {
      if (!is.null(fun)) simdata=fun(sampsize,dimension, ...)
      else simdata=rnorm(sampsize*dimension)
      nboot=nbsafe*subsampsizes[s]^dimension
      cop=corc0(simdata,sampsize,dimension,subsampsizes[s],nboot,42)
      tailcop=subsampsizes[s]^dimension
      nbootreel=cop[tailcop+2]
      cop=cop[1:tailcop]/nbootreel
      lrs[e,s,]=distance(cop,rep(1/tailcop,tailcop))
      lrs2mean[e,s,]=distance(cop,vraicop[[s]])
      scarcities[e,s]=sum(cop==0)/(nbootreel*subsampsizes[s])
    }
  }
  return(list(lrs=lrs,lrs2mean=lrs2mean,scarcities=scarcities,DistTypes=DistTypes))
}