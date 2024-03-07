#mixnorm distribution


dmixnorm = function(x, mus, sigmas, probs){
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  mapply(function(mu_, sigma_, p_) p_*dnorm(x = x, mean = mu_, sd = sigma_), mus, sigmas, probs) %>% 
    (function(yy){ 
      if(!is.array(yy) || length(dn <- dim(yy)) < 2L) sum(yy)
      else rowSums(yy)
    } )
}

pmixnorm = function(x, mus, sigmas, probs){
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  mapply(function(mu_, sigma_, p_) p_*pnorm(q = x, mean = mu_, sd = sigma_), mus, sigmas, probs) %>% 
    (function(yy){ 
      if(!is.array(yy) || length(dn <- dim(yy)) < 2L) sum(yy)
      else rowSums(yy)
    } )
} 

qmixnorm = function(p, mus, sigmas, probs, interval = c(-10**6, 10**6)){
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  sapply(p, FUN = function(pp){
    uniroot(f = function(x) {pmixnorm(x, mus, sigmas, probs)-pp}, 
            interval = interval, 
            extendInt = 'yes')$root
  })
}

normalise_weights = function(unnormalized_weights){
  unnormalized_weights / sum(unnormalized_weights)
}

rmixnorm = function(n, mus, sigmas, probs) {
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  rmultinom(n = 1, size = n, prob = probs) %>%
    as.numeric() %>%
    (function(x) {
      mapply(FUN = function(n_, mean_, sd_) {rnorm(n_, mean_, sd_)}, x, mus, sigmas)
    }) %>%
    unlist() %>% 
    as.vector()
}
# 
# 
# rmixnorm = function(n, means, sds, unnormalized_weights = NULL) {
#   if (is.null(unnormalized_weights)){
#     weights = rep(1, length(means)) / length(means)
#   }
#   else{
#     weights = unnormalized_weights  %>% normalise_weights
#   }
#   if (length(unique(length(means), length(sds), length(weights))) != 1){
#     stop('Check your number of means, sds, unnormalized_weights')
#   }
#   
#   rmultinom(n = 1, size = n, prob = weights) %>%
#     as.numeric() %>%
#     (function(x) {
#       mapply(FUN = function(n_, mean_, sd_) {rnorm(n_, mean_, sd_)}, x, means, sds)
#     }) %>%
#     unlist() %>% 
#     as.vector()
# }
