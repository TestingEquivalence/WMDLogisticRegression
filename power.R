randomExteriorPoint<-function(p,mdr, eps){
  repeat{
    skip= FALSE
    sp=resample.p(mdr$weights,p)
    mdr= updateMinDistanceModel(sp,mdr)
    
    if (mdr$min.distance>=eps){
      return(sp)
    }
  }
}

linearBoundaryPoint<-function(interiorPoint,exteriorPoint,mdr, eps){
  mdr$test=asymptotic
  aim<-function(a){
    #print(a)
    lc=a*interiorPoint+(1-a)*exteriorPoint
    nmdr=updateMinDistanceModel(lc,mdr)
    return(nmdr$min.distance-eps)
  }
  
  a=uniroot(aim, c(0,1))$root
  lc=a*interiorPoint+(1-a)*exteriorPoint
  nmdr=updateMinDistanceModel(lc,mdr)
  return(lc)
}


simulatePowerAtBoundary<-function(p,mdr, nSimulation, eps){
  set.seed(01032020)
  exteriorPoints=list()
  bndPoints=list()
  test=mdr$test
  mdr$test=asymptotic
  nPoints=100
  
  for (i in c(1:nPoints)){
    exteriorPoints[[i]]=randomExteriorPoint(p,mdr,eps)
  }
  
  for (i in c(1:nPoints)){
    bndPoints[[i]]=linearBoundaryPoint(interiorPoint = p,
                                       exteriorPoint = exteriorPoints[[i]],
                                       mdr, eps)
  }
  
  mdr$test=test

  cl=getCluster()
  power=parSapply(cl,bndPoints, simulatePowerAtPoint,mdr, nSimulation,eps)
  stopCluster(cl)

  # power=sapply(bndPoints, simulatePowerAtPoint,mdr, nSimulation,eps)
  
  return(power)
  }
  
  

