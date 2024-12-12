source("asymptoticTest.R")
source("asymptoticTestBootstrapVariance.R")
source("empiricalBootstrapTest.R")
source("simulation.R")
source("BootstrapTestTPercentile.R")
library(minpack.lm)
library(parallel)

logit=qlogis
logistic = plogis
none=""
asymptotic="asymptotic"
asymptoticBootstrapVariance="asymptoticBootstrapVariance"
empiricalBootstrap="empiricalBootstrap"
tPercentileBootstrap="tPercentileBootstrap"

min_dst_logit<-function(formula,data, weights,  test, alpha,
                        nSimulation, fixIntercept){
  
  #initial information
  mdr=list()
  mdr$data=data
  mdr$formula=as.formula(formula)
  mdr$frm=formula
  mdr$n=sum(weights)
  mdr$weights=weights
  mdr$alpha=alpha
  mdr$test=test
  mdr$nSimulation=nSimulation
 
  #logit regression for initial values
  lr <- glm(mdr$formula,mdr$data, family = quasibinomial("logit"), weights =mdr$weights)
  
  # dummy model for technical reasons
  md= lm(mdr$frm, mdr$data)
  y=all.vars(as.formula(mdr$frm))[1]
  
  # logistic model for given parameters
  w=mdr$weights/mdr$n
  
  # fix intercept if necessary
  intercept=lr$coefficients[[1]]
  
  distance<-function(par){
    md$coefficients=getCoef(par,intercept, fixIntercept)
    l=predict.lm(md,mdr$data)
    # print(par)
    (logistic(l)-mdr$data[[y]])*w
  }
  
  # calculate minimum distance estimator
  if (fixIntercept){
    par=lr$coefficients[-1]
  }
  else{
    par=lr$coefficients
  }
  res=nls.lm(par=par, fn=distance)
  mdr$result.nls.lm=res
  mdr$fixIntercept=fixIntercept
  
  # calculate min distance
  mdr$min.distance=sqrt(deviance(res))
  par=coef(res)
  mdr$coefficients=getCoef(par,intercept,fixIntercept)
 
  # calculate fitted
  md$coefficients=mdr$coefficients
  l=predict.lm(md,mdr$data)
  mdr$fitted=logistic(l)
   
  # easy access to other data 
  mdr$y=mdr$data[[y]]
  mdr$w=mdr$weights/mdr$n
  
  # test results
  mdr$min.epsilon=NA
  
  if (asymptotic==test) {
    mdr$min.epsilon=asymptoticTest(mdr=mdr)
  }
  
  if (asymptoticBootstrapVariance==test){
    mdr$min.epsilon=asymptoticTestBootstrapVariance(mdr,nSimulation)
  }
  
  if (empiricalBootstrap==test){
    mdr$min.epsilon=empiricalBootstrapTest(mdr,nSimulation)
  }
  
  if (tPercentileBootstrap==test){
    mdr$min.epsilon=tPercentileBootstrapTest(mdr,nSimulation)
  }
  
  return(mdr)
}

updateMinDistanceModel<-function(p,weights,mdr){
  df=mdr$data
  y=all.vars(as.formula(mdr$frm))[1]
  df[[y]]=p
  
  nlr=min_dst_logit(mdr$frm,data=df,weights=weights,test = mdr$test, alpha = mdr$alpha,
                    nSimulation = mdr$nSimulation, fixIntercept = mdr$fixIntercept)
  return(nlr)
}

getCoef<-function(par, intercept, fixIntercept){
   if (fixIntercept){
     npar=names(par)
     npar=c("Intercept",npar)
     par=c(intercept,par)
     names(par)=npar
     return(par)
   }
  return(par)
}
  