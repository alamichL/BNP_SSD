library(Rcpp)
sourceCpp("scripts/get_quantiles_BNP_density.cpp")

fill_sigmas = function(fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, fit$means, fit$sigmas)
}

# get_median_quantile_BNP_density = function(fit, q = 0.05, search_interval = c(-10**6, 10**6)){
#   get_median_quantile(q = q, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights, inf_limit = search_interval[[1]], sup_limit = search_interval[[2]])
# }


get_median_quantile_full_BNPdensity = function(fit, q = 0.05, search_interval = c(-10**6, 10**6)){

  get_median_quantile(q = q, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights, inf_limit = search_interval[[1]], sup_limit = search_interval[[2]])
}

get_median_quantile_semi_BNPdensity = function(fit, q = 0.05, search_interval = c(-10**6, 10**6)){

  fit$sigmas_filled = fill_sigmas(fit)

  get_median_quantile(q = q, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights, inf_limit = search_interval[[1]], sup_limit = search_interval[[2]])
}

get_median_quantile_and_CI_full_BNP_density = function(fit, q = 0.05, search_interval = c(-10**6, 10**6)){
  get_median_quantile_and_CI_BNP(q = q, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights, inf_limit = search_interval[[1]], sup_limit = search_interval[[2]])
}

get_median_quantile_and_CI_semi_BNP_density = function(fit, q = 0.05, search_interval = c(-10**6, 10**6)){

  fit$sigmas_filled = fill_sigmas(fit)

  get_median_quantile_and_CI_BNP(q = q, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights, inf_limit = search_interval[[1]], sup_limit = search_interval[[2]])
}
