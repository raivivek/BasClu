#include "RcppArmadillo.h"
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace std;

// [[Rcpp::export]]
List test(arma::mat& i) {
  return Rcpp::wrap(i.each_col());
}
