

asymptStDev<-function(mdr){
  
  w=mdr$weights/mdr$n
  
  y=all.vars(as.formula(mdr$frm))[1]
  p=mdr$data[[y]]
  
  #calculate q  vector
  q1=p*w
  q0=(1-p)*w
  q=c(q0,q1)
  
  #calculate derivative
  r=mdr$residuals
  dq0=r*r-2*r*p
  dq1=r*r+2*r*q0/(q1+q0)
  
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
