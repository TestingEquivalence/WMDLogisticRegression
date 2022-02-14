source("asymptoticTest.R")
source("asymptoticTestBootstrapVariance.R")
source("empiricalBootstrapTest.R")
source("simulation.R")
source("BootstrapTestTPercentile.R")
library(minpack.lm)

logit=qlogis
logistic = plogis
none=""
asymptotic="asymptotic"
asymptoticBootstrapVariance="asymptoticBootstrapVariance"
empiricalBootstrap="empiricalBootstrap"
tPercentileBootstrap="tPercentileBootstrap"

min_dst_logit<-function(formula,data, weights,  test, alpha=0.05,
                        nSimulation=200){
  
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
  distance<-function(coef){
    md$coefficients=coef
    l=predict.lm(md,mdr$data)
    (logistic(l)-mdr$data[[y]])*w
  }
  
  # calculate minimum distance estimator
  
  res=nls.lm(par=lr$coefficients, fn=distance)
  mdr$result.nls.lm=res
  
  # calculate min distance
  mdr$min.distance=sqrt(deviance(res))
  mdr$coefficients=coef(res)
  mdr$residuals=res$fvec
  mdr$fitted=res$fvec+mdr$data[[y]]
  
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

updateMinDistanceModel<-function(p,mdr){
  df=mdr$data
  y=all.vars(as.formula(mdr$frm))[1]
  df[[y]]=p
  
  nlr=min_dst_logit(mdr$frm,data=df,weights=mdr$weights,test = mdr$test, alpha = mdr$alpha,
                    nSimulation = mdr$nSimulation)
  return(nlr)
}
