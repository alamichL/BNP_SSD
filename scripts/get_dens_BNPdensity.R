
fill_sigmas = function(fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, fit$means, fit$sigmas)
}

library(Rcpp)




cppFunction('NumericVector dmixnorm_vec_loop(NumericVector xs, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list) {

  int nit = means_list.size(); 
  NumericVector res(xs.size());
  //std::fill(res.begin(), res.end(), 0);

  
  for(int i = 0; i < nit; ++i) {
    NumericVector mus =  as<NumericVector>(means_list[i]);
    NumericVector sigmas =  as<NumericVector>(sigmas_list[i]);
    NumericVector weights =  as<NumericVector>(weights_list[i]);

    for(int j = 0; j < mus.size(); ++j) {
      res += dnorm(xs, mus[j], sigmas[j])*weights[j];
    }
    
  }
  
  return res*1./nit;
}')


get_dens_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  
  dmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights)
}

get_dens_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  
  fit$sigmas_filled = fill_sigmas(fit)
  
  dmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights)
}  
  