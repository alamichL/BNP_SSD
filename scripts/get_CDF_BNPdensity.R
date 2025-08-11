fill_sigmas = function(fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, fit$means, fit$sigmas)
}


cppFunction('NumericVector pmixnorm_vec_loop(NumericVector xs, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list) {

  int nit = means_list.size(); 
  NumericVector res(xs.size());
  //std::fill(res.begin(), res.end(), 0);

  
  for(int i = 0; i < nit; ++i) {
    NumericVector mus =  as<NumericVector>(means_list[i]);
    NumericVector sigmas =  as<NumericVector>(sigmas_list[i]);
    NumericVector weights =  as<NumericVector>(weights_list[i]);

    for(int j = 0; j < mus.size(); ++j) {
      res += pnorm(xs, mus[j], sigmas[j])*weights[j];
    }
    
  }
  
  
  return res*1./nit;
}')


get_CDF_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  
  pmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights)
  
}

get_CDF_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  
  fit$sigmas_filled = fill_sigmas(fit)
  
  pmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights)
  
}  

cppFunction('NumericMatrix pmixnorm_matrix(NumericVector xs, 
                                           Rcpp::List means_list, 
                                           Rcpp::List sigmas_list, 
                                           Rcpp::List weights_list) {
  int nit = means_list.size(); 
  int nx = xs.size();
  NumericMatrix res(nit, nx);

  for(int i = 0; i < nit; ++i) {
    NumericVector mus = as<NumericVector>(means_list[i]);
    NumericVector sigmas = as<NumericVector>(sigmas_list[i]);
    NumericVector weights = as<NumericVector>(weights_list[i]);

    for(int j = 0; j < mus.size(); ++j) {
      for(int k = 0; k < nx; ++k) {
        res(i, k) += R::pnorm(xs[k], mus[j], sigmas[j], 1, 0) * weights[j];
      }
    }
  }

  return res;
}')

get_CDF_with_quantiles_BNPdensity <- function(fit, xs = seq(-5, 5, length.out = 100), probs = c(0.025, 0.975)) {
  # Fill sigmas if needed
  if (is.null(fit$sigmas_filled)) {
    fit$sigmas_filled <- fill_sigmas(fit)
  }
  
  # Get matrix of CDF values: rows = iterations, cols = xs
  cdf_matrix <- pmixnorm_matrix(xs = xs, 
                                means_list = fit$means, 
                                sigmas_list = fit$sigmas_filled, 
                                weights_list = fit$weights)
  
  # Compute pointwise mean and quantiles
  mean_cdf <- colMeans(cdf_matrix)
  quantiles_cdf <- apply(cdf_matrix, 2, quantile, probs = probs)
  
  # Return as a data frame
  result <- data.frame(xs = xs,
                       mean = mean_cdf,
                       lower = quantiles_cdf[1, ],
                       upper = quantiles_cdf[2, ])
  return(result)
}
