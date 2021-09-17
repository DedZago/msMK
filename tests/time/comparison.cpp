#include <RcppArmadillo.h>
#include "bintree.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// Tony
// [[Rcpp::export]]
void bintree_tony(int smax){
  int nelem = (int) std::pow(2.0, smax + 1) - 1;
  NumericVector v(nelem, 0.0);
  bintree* root = new struct bintree;
  setTree(0.0, root);
  array2tree(REAL(v), smax, root);
  delete root;
};

// Dede 
// [[Rcpp::export]]
void bintree_dede(int smax){
  int nelem = (int) std::pow(2.0, smax + 1) - 1;
  field<double> vec(nelem);
  for(int i = 0; i < nelem; ++i){
    vec(i) = 0.0;
  }
};


// [[Rcpp::export]]
double rtmvnorm_arma_botev(const double& mean,
                           const double& sigma,
                           const double& lower,
                           const double& upper){
  Environment env("package:TruncatedNormal"); 
  Function fun = env["rtmvnorm"];
  
  NumericMatrix Sigma(1, 1);
  Sigma.fill(sigma);
  NumericVector out = fun(Rcpp::_["n"]  = 1,
                          Rcpp::_["mu"] = mean,
                          Rcpp::_["sigma"] = Sigma,
                          Rcpp::_["lb"] = lower, 
                          Rcpp::_["ub"] = upper); 
  return out(0);
}

// [[Rcpp::export]]
void lol(int n){
  double a = 0.0;
  for(int i = 0; i < n; ++i){
    a = a + 1;
  }
}

/*** R

rtmvnorm_arma_botev(0, 1, 0, 5)

*/
