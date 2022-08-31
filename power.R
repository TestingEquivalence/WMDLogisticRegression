randomExteriorPoint<-function(p,mdr, eps){
  repeat{
    skip= FALSE
    #resample cell sizes first
    nn=rmultinom(1, mdr$n,mdr$weights)
    #resample counting frequencies
    np=resample.p(nn,p)
    mdr= updateMinDistanceModel(np,nn,mdr)
    
    if (mdr$min.distance>=eps){
      return(mdr)
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
    new_p=new_count/new_weights
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


simulatePowerAtBoundary<-function(p,mdr, nSimulation, eps){
  set.seed(01032020)
  exteriorModels=list()
  bndModels=list()
  test=mdr$test
  mdr$test=none
  nPoints=100
  
  for (i in c(1:(2*nPoints))){
    exteriorModels[[i]]=randomExteriorPoint(p,mdr,eps)
  }
  
  j=1
  i=1
  while(i<=nPoints){
    tryCatch({
      bndModels[[i]]=linearBoundaryPoint(mdr,exteriorModels[[j]],eps)
      i=i+1
    }, 
    error=function(e){
      print("error")
    },finally = {
    })
    print(j)
    j=j+1
  }
  
  mdr$test=test

  # cl=getCluster()
  # power=parSapply(cl,bndPoints, simulatePowerAtPoint,mdr, nSimulation,eps)
  # stopCluster(cl)

  power=sapply(bndModels, simulatePowerAtPoint, nSimulation,eps)
  
  return(power)
}
  
  

