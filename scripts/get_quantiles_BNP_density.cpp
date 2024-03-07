
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>
using namespace Rcpp;
#include <iostream>
#include <sstream>
#include <string>

NumericVector quantileCpp(NumericVector x, NumericVector probs) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int npr = probs.size();
  NumericVector ans(npr);
  for(int i=0; i<npr; i++){
    ans[i] = as<double>(quantile(x, probs[i]));
  }
return ans;
}

double median_rcpp(NumericVector x) {
   NumericVector y = clone(x);
   int n, half;
   double y1, y2;
   n = y.size();
   half = n / 2;
   if(n % 2 == 1) {
      // median for odd length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      return y[half];
   } else {
      // median for even length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      y1 = y[half];
      std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
      y2 = y[half-1];
      return (y1 + y2) / 2.0;
   }
}

NumericVector pmixnorm(NumericVector x, NumericVector mus, NumericVector sigmas, NumericVector weights) {

  NumericVector res(x.size());
  std::fill(res.begin(), res.end(), 0);
    //std::cout << "x " << x << std::endl;


    for(int j = 0; j < mus.size(); ++j) {
    //std::cout << "res " << res << std::endl;

      res += pnorm(x, mus[j], sigmas[j], true, false) * weights[j];
      //res += pnorm(x, mus[j], sigmas[j]) * weights[j];
    }
  
  return res;
}

double pmixnorm2(double x, NumericVector mus, NumericVector sigmas, NumericVector weights) {

  double res = 0;
    //std::cout << "x " << x << std::endl;


    for(int j = 0; j < mus.size(); ++j) {
    //std::cout << "res " << res << std::endl;

      res += R::pnorm(x, mus[j], sigmas[j], true, false) * weights[j];
      //res += pnorm(x, mus[j], sigmas[j]) * weights[j];
    }
  
  return res;
}


class Function_to_root {
public:
    Function_to_root(double q, NumericVector mus, NumericVector sigmas, NumericVector weights) : q(q), mus(mus), sigmas(sigmas), weights(weights) {}

    double operator()(const double x) const {
    //std::cout << "q " << q << std::endl;

        return pmixnorm2(x, mus, sigmas, weights)  - q;
    }
private:
    double q;
    NumericVector mus;
    NumericVector sigmas;
    NumericVector weights;
};

double find_the_quantile(double q, NumericVector mus, NumericVector sigmas, NumericVector weights, double inf_limit, double sup_limit) {
    //std::cout << "q " << q << std::endl;

    Function_to_root t( q, mus, sigmas, weights);
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);

    Result r1 = boost::math::tools::toms748_solve(t, inf_limit, sup_limit, tol, max_iter);
    //std::cout << "root bracketed: [ " << r1.first << " , " << r1.second <<  " ]" << std::endl;
    //std::cout << "f("<< r1.first << ")=" << t(r1.first) << std::endl;
    //std::cout << "f("<< r1.second << ")=" << t(r1.second) << std::endl;
    //std::cout << "max_iter=" << max_iter << std::endl;
    return r1.first;
}

NumericVector get_all_quantiles(double q, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list, double inf_limit, double sup_limit) {
  //std::cout << "q " << q << std::endl;
  
  // std::cout << "inf_limit=" << inf_limit << "sup_limit=" << sup_limit << std::endl;
  
  
  int nit = means_list.size(); 
  NumericVector res(nit);
  
  for(int i = 0; i < nit; ++i) {
    NumericVector mus =  as<NumericVector>(means_list[i]);
    NumericVector sigmas =  as<NumericVector>(sigmas_list[i]);
    NumericVector weights =  as<NumericVector>(weights_list[i]);
    res[i] = find_the_quantile(q, mus, sigmas, weights, inf_limit, sup_limit);
    //std::cout << "res=" << res << std::endl;
    
  }
  
  return res;
}


// [[Rcpp::export]]
double get_median_quantile(double q, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list, double inf_limit, double sup_limit) {

    
  return median_rcpp(get_all_quantiles(q, means_list, sigmas_list, weights_list, inf_limit, sup_limit));
}

// [[Rcpp::export]]
NumericVector get_median_quantile_and_CI_BNP(double q, Rcpp::List means_list, Rcpp::List sigmas_list, Rcpp::List  weights_list, double inf_limit, double sup_limit) {
    
    NumericVector probs = NumericVector::create(0.5, 0.025, 0.975);
  
    return quantileCpp(get_all_quantiles(q, means_list, sigmas_list, weights_list, inf_limit, sup_limit), probs);
}

