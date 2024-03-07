#emulation of ETX

etx=function(dat,hcx=0.05){
  ###takes data and hcx as input, returns the hcx and its 90% confidence interval
  ###Fully self-contained version
  ldat=log10(dat) ####log transformation of the data
  
  ndat=length(ldat)
  
  ##### moment calculation for moment matching
  mu=mean(ldat) 
  sig=sqrt(var(ldat))
  
  
  #####function which has to be inversed to calculate the interpolation constant
  
  fzksnct=function(ks,n,conf,kp){
    x=ks*sqrt(n)
    c=kp*sqrt(n)
    pt(q=x,df=n-1,ncp=c)-conf   ######student non central distribution
  }
  
  
  kp=qnorm(hcx)
  
  ####calculation of the constants
  ksmin = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.05,kp=kp)$root
  ksmed = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.5,kp=kp)$root
  ksmax = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.95,kp=kp)$root
  
  
  
  #result
  list('HC'=10^(mu+sig*ksmed),'C_inf'=10^(mu+sig*ksmin),'C_sup'=10^(mu+sig*ksmax))
}


etx_on_log=function(ldat,hcx=0.05){
  
  ###takes data and hcx as input, returns the hcx and its 90% confidence interval
  
  ndat=length(ldat)
  
  ##### moment calculation for moment matching
  mu=mean(ldat) 
  sig=sqrt(var(ldat))
  
  
  
  #####function which has to be inversed to calculate the interpolation constant
  
  fzksnct=function(ks,n,conf,kp){
    x=ks*sqrt(n)
    c=kp*sqrt(n)
    pt(q=x,df=n-1,ncp=c)-conf   ######student non central distribution
  }
  
  
  kp=qnorm(hcx)
  
  ####calculation of the constants
  ksmin = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.05,kp=kp)$root
  ksmed = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.5,kp=kp)$root
  ksmax = uniroot(fzksnct,interval=c(-10,100),n=ndat,conf=0.95,kp=kp)$root
  
  #result
  list('HC'=mu+sig*ksmed,'C_inf'=mu+sig*ksmin,'C_sup'=mu+sig*ksmax)
}

# etx_on_log_fast--------------------------------------------------------------------------
# Precomputation

library(dplyr)
library(parallel)

ndatmax = 1000


#####function which has to be inversed to calculate the interpolation constant

fzksnct=function(ks,n,conf,kp){
  x=ks*sqrt(n)
  c=kp*sqrt(n)
  pt(q=x,df=n-1,ncp=c)-conf   ######student non central distribution
}


kp=qnorm(0.05)

####calculation of the constants
ksmin = 3:ndatmax %>% 
  mclapply(FUN = function(i) uniroot(fzksnct,interval=c(-10,100),n=i,conf=0.05,kp=kp)$root, mc.cores = detectCores()) %>% 
  unlist()
ksmed = 3:ndatmax %>% 
  mclapply(FUN = function(i) uniroot(fzksnct,interval=c(-10,100),n=i,conf=0.5,kp=kp)$root, mc.cores = detectCores()) %>% 
  unlist()
ksmax = 3:ndatmax %>% 
  mclapply(FUN = function(i) uniroot(fzksnct,interval=c(-10,100),n=i,conf=0.95,kp=kp)$root, mc.cores = detectCores()) %>% 
  unlist()

etx_on_log_fast=function(ldat,hcx=0.05){
  
  
  ###takes data and hcx as input, returns the hcx and its 90% confidence interval
  
  ndat=length(ldat)
  
  
  if(hcx!=0.05|ndat>ndatmax) etx_on_log(ldat,hcx=0.05)
  
  else{
  ##### moment calculation for moment matching
  mu=mean(ldat) 
  sig=sqrt(var(ldat))
  
  
  #result
  list('HC'=mu+sig*ksmed[[ndat]],'C_inf'=mu+sig*ksmin[[ndat]],'C_sup'=mu+sig*ksmax[[ndat]])
  }
}