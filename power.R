randomExteriorPoint<-function(p,mdr, eps){
  repeat{
    #resample cell sizes first
    nn=rmultinom(1, mdr$n,mdr$weights)
    #resample counting frequencies
    np=resample.p(nn,p)
    nmdr= updateMinDistanceModel(np,nn,mdr)
    # print(nmdr$min.distance)
    
    if (nmdr$min.distance>=eps){
      return(nmdr)
    }
  }
}

linearBoundaryPoint<-function(int_mdr,ext_mdr, eps){
  mdr$test=asymptotic
  linear_mdr<-function(a){
    #print(a)
    new_weights=a*int_mdr$weights+(1-a)*ext_mdr$weights
    int_count=int_mdr$y*int_mdr$weights
    ext_count=ext_mdr$y*ext_mdr$weights
    new_count=a*int_count+(1-a)*ext_count
    new_p= ifelse(new_weights==0,0, new_count/new_weights)
    nmdr=updateMinDistanceModel(new_p,new_weights,int_mdr)
    return(nmdr)
  }
  
  aim<-function(a){
    nmdr=linear_mdr(a)
    return(nmdr$min.distance-eps)
  }
  
  a=uniroot(aim, c(0,1))$root
  nmdr=linear_mdr(a)
  return(nmdr)
}

generateBoundaryPoints<-function(p,mdr,nSimulation, eps){
  set.seed(01032020)
  exteriorModels=list()
  bndModels=list()
  nPoints=100
  
  # set test to none to increase performance 
  # by searching for the boundary points
  test=mdr$test
  mdr$test=none
  
  
  for (i in c(1:(nPoints+10))){
    exteriorModels[[i]]=randomExteriorPoint(p,mdr,eps)
    print(i)
  }
  
  j=1
  i=1
  while(i<=nPoints){
    tryCatch({
      bmdr=linearBoundaryPoint(mdr,exteriorModels[[j]],eps)
      # set test back to the correct value
      bmdr$test=test
      
      param=list()
      param$mdr=bmdr
      param$nSimulation=nSimulation
      param$eps=eps
      param$nr=i
      
      bndModels[[i]]=param
      print("ok")
      print(i)
      i=i+1
    }, 
    error=function(e){
      print("error")
    },finally = {
    })
    
    print(j)
    j=j+1
  }
  
  return(bndModels)
}


simulatePowerAtBoundary<-function(bndPoints){
  set.seed(12022025)
  cl=getCluster()
  power=parSapply(cl,bndPoints, simulatePowerAtPoint)
  stopCluster(cl)

  # power=sapply(bndModels, simulatePowerAtPoint)
  
  # power=rep(0,nPoints)
  # for (i in c(1:nPoints)){
  #   power[i]=simulatePowerAtPoint(bndModels[[i]], nSimulation,eps)
  #   print(i)
  # }
  
  nPoints=length(bndPoints)
  for (i in c(1:nPoints)){
    fname=paste0("r",i,".csv")
    file.remove(fname)
  }
  return(power)
}
  
  

