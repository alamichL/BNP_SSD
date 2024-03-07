#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
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
NumericVector rowMedianRcpp(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  int position = nrow / 2; // Euclidian division
  NumericVector out(nrow);
  for (int j = 0; j < nrow; j++) { 
    NumericVector y = x(j,_); // Copy row -- original will not be mod
    std::nth_element(y.begin(), y.begin() + position, y.end()); 
    out[j] = y[position];  
  }
  return out;
}