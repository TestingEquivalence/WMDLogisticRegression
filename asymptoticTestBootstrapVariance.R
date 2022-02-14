bootstrapVolatility<-function(mdr,nSimulation){
  #calculate bootstrap volatility
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
 
  mdr$test=none
  res=rep(NA,nSimulation)
  
  for (i in c(1:nSimulation)){
    #resample cell sizes first
    n=rmultinom(1, mdr$n,mdr$weights)
    np=resample.p(n,p)
    nmdr=updateMinDistanceModel(p=np,mdr)
    res[i]=nmdr$min.distance^2
  }

    return(sd(res))
}

asymptoticTestBootstrapVariance<-function(mdr,nSimulation){
  #calculate asymptotic min eps
  vol = bootstrapVolatility(mdr,nSimulation)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
