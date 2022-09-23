#source("asymptoticTest.R")
source("simulation.R")

simulatePowerMDR<-function(p, nSimulation, mdr){
  set.seed(01032020)
  nsim=list()
  psim=list()
  for (i in c(1:nSimulation)){
    #resample cell sizes first
    nsim[[i]]=rmultinom(1, mdr$n,mdr$weights)
    #resample counting frequencies
    psim[[i]]=resample.p(nsim[[i]],p)
  }
  
  res=list()
  for (i in c(1:nSimulation)){
    nlr=updateMinDistanceModel(psim[[i]],nsim[[i]],mdr)
    res[[i]]=mdr2results(nlr)
  }

  return(res)
}

simulatePowerLR<-function(p, nSimulation, lr, mdr){
  set.seed(01032020)
  nsim=list()
  psim=list()
  for (i in c(1:nSimulation)){
    #resample cell sizes first
    nsim[[i]]=rmultinom(1, mdr$n,mdr$weights)
    #resample counting frequencies
    psim[[i]]=resample.p(nsim[[i]],p)
  }
  
  res=list()
  for (i in c(1:nSimulation)){
    nlr=updateLogitModel(psim[[i]],nsim[[i]],lr)
    res[[i]]=mdr2results(nlr)
  }
  
  return(res)
}
