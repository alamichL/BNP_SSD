
dt_ = function(x, df, mean, sd){
  dt((x-mean)/sd, df, ncp = 0)/sd
}

pt_ = function(x, df, mean, sd){
  pt((x-mean)/sd, df, ncp = 0)
}

qt_ = function(p, df, mean, sd){
  sd*qt(p, df, ncp = 0) + mean
}

rt_ = function(n, df, mean, sd){
  mean + sd*rt(n, df, ncp = 0)
}

dhalft = cmpfun(function(x, df = 1, mean = 0, sd = 1){
  dt_(x, df, mean, sd)/(1-pt_(0, df, mean, sd))
})

phalft = cmpfun(function(q, df = 1, mean = 0, sd = 1){
  (pt_(q, df, mean, sd)-pt_(0, df, mean, sd))/(1-pt_(0, df, mean, sd))
})

qhalft = cmpfun(function(p,  df = 1, mean = 0, sd = 1){
  qt_(p*(1-pt_(0, df, mean, sd))+pt_(0, df, mean, sd), df, mean, sd)
})

rhalft = cmpfun(function(n,  df = 1, mean = 0, sd = 1){
  abs(rt_(n,  df, mean, sd))
})