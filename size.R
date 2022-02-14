#source("asymptoticTest.R")
source("simulation.R")

simulatePowerAtModel<-function(df,n,p,lr, updateLR, nSimulation){
  set.seed(01032020)
  psim=list()
  for (i in c(1:nSimulation)){
    psim[[i]]=resample.p(n,p)
  }
  
  res=list()
  for (i in c(1:nSimulation)){
    nlr=updateLR(p=psim[[i]],lr)
    res[[i]]=mdr2results(nlr)
  }

  return(res)
}
