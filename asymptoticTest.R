

asymptStDev<-function(mdr){
  
  w=mdr$w
  p=mdr$y
  f=mdr$fitted
  
  #calculate q  vector
  q1=p*w
  q0=(1-p)*w
  
  #calculate derivative
  r=(p-f)*w
  dq0=-2*f*r
  dq1=2*(1-f)*r
  
  vol=asymptSDMultinomial(p=c(q0,q1), derivative=c(dq0,dq1))
  return(vol)
}

asymptSDMultinomial<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  
  f<-function(j){
    v=vec[j]*vec
    v=v*p[j]
    v=v*p
    return(sum(v))
  }
  
  vv=sapply(c(1:k),f)
  vnsq_2=sum(vv)
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}



asymptoticTest<-function(mdr){
  #calculate asymptotic min eps
  vol = asymptStDev(mdr)/sqrt(mdr$n)
  qt=qnorm(1-mdr$alpha,0,1)
  aps = mdr$min.distance^2 + qt*vol
  return(sqrt(aps))
}
