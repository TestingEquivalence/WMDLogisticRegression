library(parallel)

resample.p<-function(n,p){
  f<-function(r){
    if (r[1]>0){
      return(rbinom(1,r[1],r[2])/r[1])
    }
    else{
      return(0)
    }
  }
  
  pdf=data.frame(n=n,p=p)
  res=apply(X=pdf, MARGIN=1,FUN = f)
  res=as.vector(res)
  return(res)
}

mdr2results<-function(mdr){
  #gather resuts for effcience purpose
  res=list(min.distance=mdr$min.distance,
           min.epsilon=mdr$min.epsilon,
           coefficients=mdr$coefficients
           )

}

getTitle<-function(res){
  #title
  title=paste("min_distance"," min_epsilon",sep=",")
  
  for (cf in names(res[[3]])) {
    title=paste(title,cf,sep=",")
  }
  
  return(title)
}

res2row<-function(row){
  sl=paste(row[[1]],row[[2]], sep=",")
  for (cf in row[[3]]){
    sl=paste(sl,cf,sep=",")
  }
  return(sl)
}

write.result<-function(mdr,fname){
  res=mdr2results(mdr)
  fc=file(fname)
  title=getTitle(res)
  row=res2row(res)
  rows=c(title,row)
  writeLines(rows,con=fc)
  close(fc)
}

write.results<-function(res,fname){
    fc=file(fname)
    title=getTitle(res[[1]])
    
    #data
    rows=sapply(res, res2row)
    rows=c(title,rows)
    #rows=as.vector(rows)
    
    writeLines(rows,con=fc)
    close(fc)
}

updateLogitModel<-function(p,lr){
  df=lr$data
  frm=lr$formula
  depVar=all.vars(as.formula(frm))[1]
  
  df[[depVar]]=p
  nlr=update(lr,data=df)
  v=(p-nlr$fitted.values)*weights(lr)/sum(weights(lr))
  nlr$min.distance=sqrt(sum((v)^2))
  if (is.infinite(nlr$min.distance)){browser()}
  
  return(nlr)
}


simulatePowerAtPoint<-function(p,mdr, nSimulation,eps){
  set.seed(01032020)
  psim=list()
  for (i in c(1:nSimulation)){
    psim[[i]]=resample.p(mdr$weights,p)
  }
  
  res=rep(NA,nSimulation)
  
  for (i in c(1:nSimulation)){
    nmdr=updateMinDistanceModel(p=psim[[i]],mdr)
    res[i]=nmdr$min.epsilon
  }
  
  return(sum(res<=eps)/nSimulation)
}


# Calculate the number of cores
getCluster<-function(){
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores,'SOCK')
  clusterExport(cl,c("min_dst_logit","resample.p","updateMinDistanceModel","simulatePowerAtPoint",
                     "logit", "logistic","asymptStDev","asymptoticTest","none",
                     "linearBoundaryPoint","nls.lm","asymptotic","asymptSDMultinomial",
                     "asymptoticBootstrapVariance","empiricalBootstrap","bootstrapVolatility",
                     "asymptoticTestBootstrapVariance","empiricalBootstrapTest",
                     "tPercentileBootstrapTest","tPercentileBootstrap"))
  
  return(cl)
}