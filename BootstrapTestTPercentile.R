tPercentileBootstrapTest<-function(mdr,nSimulation){
  #calculate bootstrap distribution
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
  mdr$test=none
  
  res=rep(NA,nSimulation)
  stDev=asymptStDev(mdr)
  
  for (i in c(1:nSimulation)){
    #resample cell sizes first
    n=rmultinom(1, mdr$n,mdr$weights)
    #resample counting frequencies
    np=resample.p(n,p)
    nmdr=updateMinDistanceModel(np,n,mdr)
    nStDev=asymptStDev(nmdr)
    res[i]=(nmdr$min.distance^2-mdr$min.distance^2)/nStDev
  }
  
  #calculate quantile of bootstrap distribution
  qt=quantile(res,mdr$alpha,type=1)
  min_eps=mdr$min.distance^2-stDev*qt
  if (min_eps<0) {return(0)}
  return(sqrt(min_eps))
}