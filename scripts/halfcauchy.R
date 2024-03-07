
dhalfcauchy = cmpfun(function(x, location = 0, scale = 1){
  dcauchy(x, location, scale)/(1-pcauchy(0, location, scale))
})

phalfcauchy = cmpfun(function(q, location = 0, scale = 1){
  (pcauchy(q, location, scale)-pcauchy(0, location, scale))/(1-pcauchy(0, location, scale))
})



qhalfcauchy = cmpfun(function(p, location = 0, scale = 1){
  qcauchy(p*(1-pcauchy(0, location, scale))+pcauchy(0, location, scale), location, scale)
})


rhalfcauchy = cmpfun(function(n, location = 0, scale = 1){
  abs(rcauchy(n, location, scale))
})