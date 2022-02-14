empiricalBootstrapTest<-function(mdr,nSimulation){
  #calculate bootstrap distribution
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
  n=mdr$weights
  
  mdr$test=none
  res=rep(NA,nSimulation)
  
  for (i in c(1:nSimulation)){
    np=resample.p(n,p)
    nmdr=updateMinDistanceModel(p=np,mdr)
    res[i]=nmdr$min.distance^2
  }
  
  #calculate quantile of bootstrap distribution
  qt=quantile(res,mdr$alpha,type=1)
  min_eps=-qt+2*mdr$min.distance^2
  if (min_eps<0) {return(0)}
  return(sqrt(min_eps))
}