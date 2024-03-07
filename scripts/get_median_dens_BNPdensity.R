
fill_sigmas = function(fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, fit$means, fit$sigmas)
}

library(Rcpp)

sourceCpp('scripts/get_median_dens_BNP_density.cpp')

get_median_dens_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  dmixnorm_vec_loop_median(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights)
}

get_median_dens_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  fit$sigmas_filled = fill_sigmas(fit)

  dmixnorm_vec_loop_median(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights)
}
