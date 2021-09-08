#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "bintree.h"
#include "auxGibbs.h"
#include "gibbs.h"

using namespace Rcpp;
using namespace arma;

//' Compute msMK density on a set of points
//'
//' @param x: Matrix of size n x p, containing n points onto which to evaluate the density
//' @param TH: Field of vectors, sequence of location parameters for the density 
//' @param SIG: Field of matrices, sequence of scale parameters for the density
//' @param prob: vector of doubles, contains probabilities of each cluster
//'
//' @return two-column matrix where every row (s_i, h_i) is the cluster allocation of the i-th data point
// [[Rcpp::export]]
Rcpp::NumericVector dmsMK(
    const arma::mat &x, 
    const arma::mat &TH,
    const arma::cube &SIG,
    const arma::cube &thrs,
    const std::vector<double> &prob,
    bool indep)
{
    int smax = n_scales_tree(TH);
    
    vec out(x.n_rows, fill::zeros);
    
    for(int s = 0; s <= smax; s++){
      R_CheckUserInterrupt();
      int hmax = n_elem_scale(s);
      for (int h = 1; h <= hmax; h++) {
        vec theta = extractNode(TH, s, h);
        mat sigma = extractNode(SIG, s, h);
        vec dens = dmvnorm(x, theta, sigma);
        
        out += extractNode(prob, s, h) * dens; 
      }
    }
    return NumericVector(out.begin(), out.end());
}


