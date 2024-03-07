#include <Rcpp.h>
using namespace Rcpp;

NumericVector colMedianRcpp(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  int position = nrow / 2; // Euclidian division
  NumericVector out(ncol);
  for (int j = 0; j < ncol; j++) { 
    NumericVector y = x(_,j); // Copy column -- original will not be mod
    std::nth_element(y.begin(), y.begin() + position, y.end()); 
    out[j] = y[position];  
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dmixnorm_vec_loop_median(NumericVector xs, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list) {
  
  int nit = means_list.size(); 
  NumericVector res(xs.size());
  Rcpp::NumericMatrix densmat(nit,xs.size());
  //std::fill(res.begin(), res.end(), 0);
  
  
  for(int i = 0; i < nit; ++i) {
    NumericVector mus =  as<NumericVector>(means_list[i]);
    NumericVector sigmas =  as<NumericVector>(sigmas_list[i]);
    NumericVector weights =  as<NumericVector>(weights_list[i]);
    std::fill(res.begin(), res.end(), 0);
    
    for(int j = 0; j < mus.size(); ++j) {
      res += dnorm(xs, mus[j], sigmas[j])*weights[j];
    }
    
    densmat(i,_) = res;
    
    // Rcout << "densmat " << densmat << std::endl;
    // Rcout << "res " << res << std::endl;
  }
  
  return colMedianRcpp(densmat);
}