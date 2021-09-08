#include <math.h>
#include <R.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include "bintree.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef AUXGIBBS_H_ 
#define AUXGIBBS_H_

// int sampleC(int k, arma::vec prob);

inline int ceil(int a, int b) {
  // Ceiling function
  return (a / b + (a % b != 0));
}


arma::vec rinvgamma(R_xlen_t n, double shape, double scale);

double rinvgamma(double shape, double scale);

arma::vec rtmvnorm_arma(const arma::vec& mean,
                        const arma::mat& sigma,
                        const arma::vec& lower,
                        const arma::vec& upper,
                        const int& burnin,
                        const int& thinning); 

arma::mat rtmvnorm_arma(const int&n,
                        const arma::vec& mean,
                        const arma::mat& sigma,
                        const arma::vec& lower,
                        const arma::vec& upper,
                        const int& burnin,
                        const int& thinning); 

NumericMatrix rtmvnorm_rcpp(const int& n,
                       const NumericVector& mean,
                       const NumericMatrix& sigma,
                       const NumericVector& lower,
                       const NumericVector& upper,
                       const int& burnin,
                       const int& thinning);

arma::vec rtmvnorm_arma_botev(const arma::vec& mean,
                              const arma::mat& sigma,
                              const arma::vec& lower,
                              const arma::vec& upper);

arma::mat rtmvnorm_arma_botev(const int& n,
                              const arma::vec& mean,
                              const arma::mat& sigma,
                              const arma::vec& lower,
                              const arma::vec& upper);
    
double rtmvnorm_arma_botev(const double& mean,
                           const double& sigma,
                           const double& lower,
                           const double& upper);

double dtmvnorm_botev(const rowvec &x,
                           const arma::vec& mean,
                           const arma::mat& sigma,
                           const arma::vec& lower,
                           const arma::vec& upper);

NumericVector dtmvnorm_botev(const arma::mat &x,
                             const arma::vec& mean,
                             const arma::mat& sigma,
                             const arma::vec& lower,
                             const arma::vec& upper);

NumericVector dtmvnorm_botev(const arma::mat &x,
                             const arma::vec& mean,
                             const arma::mat& sigma,
                             const arma::vec& lower,
                             const arma::vec& upper,
                             bool indep);

arma::vec scaleProb(const vector<double> &prob);

void srTrees(vector<double> &Ssh,
             vector<double> &Rsh,
             const double &a,
             const double &b,
             const double &delta);
    
vector<double> computeProb(vector<double> Stree,
                           vector<double> Rtree,
                           const double &a,
                           const double &b);

void nrvTrees(vector<int> &n,
              vector<int> &r,
              vector<int> &v,
              const arma::umat &cluster);

void msMK_prior(arma::mat &TH,
                arma::cube &SIG,
                const arma::cube &thrs,
                const arma::colvec& mu0,
                const arma::colvec& k0,
                const arma::mat& sig0,
                const double& a,
                const double& b,
                const double& delta,
                const bool& indep);

arma::vec rowSum(const arma::mat &x);

#endif /* AUXGIBBS_H_ */