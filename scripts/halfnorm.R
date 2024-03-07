
dhalfnorm = cmpfun(function(x, mean = 0, sd = 1){
  dnorm(x, mean, sd)/(1-pnorm(0, mean, sd))
})

phalfnorm = cmpfun(function(q, mean = 0, sd = 1){
  (pnorm(q, mean, sd)-pnorm(0, mean, sd))/(1-pnorm(0, mean, sd))
})

qhalfnorm = cmpfun(function(p,  mean = 0, sd = 1){
  qnorm(p*(1-pnorm(0, mean, sd))+pnorm(0, mean, sd), mean, sd)
})


rhalfnorm = cmpfun(function(n,  mean = 0, sd = 1){
  abs(rnorm(n,  mean, sd))
})