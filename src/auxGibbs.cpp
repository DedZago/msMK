/*--------------------------------------------------------------------------
 Multiscale Mixtures of Kernels [msMK]
 auxGibbs.cpp - Auxiliary C++ functions to be used in the Gibbs sampling routine
 
 Version 0.1 of September 2021 
 2013 - Antonio Canale (antonio.canale@unipd.it)
 2021 - Daniele Zago (daniele.zago.1@studenti.unipd.it)
 --------------------------------------------------------------------------*/
#include <R.h>
#include "bintree.h"
#include "auxGibbs.h"
#include <Rmath.h>
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppDist.h>

// [[Rcpp::depends(truncnorm)]]

using namespace arma;
using namespace Rcpp;
using namespace std;


//-------------------------------------------------------------------------
// Simulate inverse gamma distribution
arma::vec rinvgamma(R_xlen_t n,double shape, double scale) {
        NumericVector temp = 1.0/Rcpp::rgamma(n, shape, 1/scale);
        return vec(temp.begin(), n, true, false);
}

//-------------------------------------------------------------------------
// Simulate inverse gamma distribution
// [[Rcpp::export]]
double rinvgamma(double shape, double scale) {
  return 1.0 / R::rgamma(shape, 1/scale);
}

//-------------------------------------------------------------------------
// Rcpp / RcppArmadillo interface for truncated multivariate normal simulation
// Original fortran function was modified in the last loop in order to account
// for column-major and row-major ordering difference between fortran and c++
// arrays.
extern "C" {
    void rtmvnormgibbscov_(const int *n, const int *d, const double *mean, const double *sigma,
                           const double *lower, const double *upper, const double *x0,
                           const int *burnin, const int *thinning, double *X);
}


NumericMatrix rtmvnorm_rcpp(const int& n,
                       const NumericVector& mean,
                       const NumericMatrix& sigma,
                       const NumericVector& lower,
                       const NumericVector& upper,
                       const int& burnin,
                       const int& thinning) {
    const int d = sigma.ncol();
    const NumericVector x0(d);
    NumericMatrix X(n, d);
    rtmvnormgibbscov_(&n, &d, REAL(mean), REAL(sigma), REAL(lower), REAL(upper), REAL(x0),
                      &burnin, &thinning, REAL(X));

    return X;
}

arma::vec rtmvnorm_arma(const arma::vec& mean,
                        const arma::mat& sigma,
                        const arma::vec& lower,
                        const arma::vec& upper,
                        const int& burnin,
                        const int& thinning) {
    int n = 1;
    const int d = sigma.n_cols;
    vec x0(d, fill::zeros);
    vec X(d, fill::zeros);
    rtmvnormgibbscov_(&n, &d, mean.memptr(), sigma.memptr(), lower.memptr(),
                      upper.memptr(), x0.memptr(), &burnin, &thinning, X.memptr());
    return X;
}

arma::mat rtmvnorm_arma(const int&n,
                        const arma::vec& mean,
                        const arma::mat& sigma,
                        const arma::vec& lower,
                        const arma::vec& upper,
                        const int& burnin,
                        const int& thinning) {
    const int d = sigma.n_cols;
    vec x0(d, fill::zeros);
    mat X(n, d, fill::zeros);
    rtmvnormgibbscov_(&n, &d, mean.memptr(), sigma.memptr(), lower.memptr(),
                      upper.memptr(), x0.memptr(), &burnin, &thinning, X.memptr());
    return X;
}

arma::vec rtmvnorm_arma_botev(const arma::vec& mean,
                              const arma::mat& sigma,
                              const arma::vec& lower,
                              const arma::vec& upper){
  Environment env("package:TruncatedNormal"); 
  Function fun = env["rtmvnorm"];    
  
  NumericVector out = fun(Rcpp::_["n"]  = 1,
                          Rcpp::_["mu"] = NumericVector(mean.begin(), mean.end()),
                          Rcpp::_["sigma"] = NumericMatrix(sigma.n_rows, sigma.n_cols, sigma.begin()),
                          Rcpp::_["lb"] = NumericVector(lower.begin(), lower.end()),
                          Rcpp::_["ub"] = NumericVector(upper.begin(), upper.end())); 
  
  
  return vec(out.begin(), out.length());
}

arma::mat rtmvnorm_arma_botev(const int& n,
                              const arma::vec& mean,
                              const arma::mat& sigma,
                              const arma::vec& lower,
                              const arma::vec& upper){
  Environment env("package:TruncatedNormal"); 
  Function fun = env["rtmvnorm"];    
  int p = mean.n_elem;
  NumericMatrix out = fun(Rcpp::_["n"]  = n,
                          Rcpp::_["mu"] = NumericVector(mean.begin(), mean.end()),
                          Rcpp::_["sigma"] = NumericMatrix(sigma.n_rows, sigma.n_cols, sigma.begin()),
                          Rcpp::_["lb"] = NumericVector(lower.begin(), lower.end()),
                          Rcpp::_["ub"] = NumericVector(upper.begin(), upper.end())); 
  
  
  return mat(out.begin(), n, p);
}

double rtmvnorm_arma_botev(const double& mean,
                           const double& sigma,
                           const double& lower,
                           const double& upper){
  Environment env("package:truncnorm"); 
  Function fun = env["rtruncnorm"];    
  
  double out = *REAL(fun(Rcpp::_["n"]  = 1,
                   Rcpp::_["mean"] = mean,
                   Rcpp::_["sd"] = sqrt(sigma),
                   Rcpp::_["a"] = lower, 
                   Rcpp::_["b"] = upper)); 
  
  return out;
}


double dtmvnorm_botev(const rowvec &x,
                      const arma::vec& mean,
                      const arma::mat& sigma,
                      const arma::vec& lower,
                      const arma::vec& upper){
  Environment env("package:TruncatedNormal"); 
  
  Function fun = env["dtmvnorm"];    
  
  NumericVector out = fun(Rcpp::_["x"]  = NumericMatrix(1, x.n_elem, x.begin()),
                          Rcpp::_["mu"] = NumericVector(mean.begin(), mean.end()),
                          Rcpp::_["sigma"] = NumericMatrix(sigma.n_rows, sigma.n_cols, sigma.begin()),
                          Rcpp::_["lb"] = NumericVector(lower.begin(), lower.end()),
                          Rcpp::_["ub"] = NumericVector(upper.begin(), upper.end())
                          ); 
  
  
  return out[0];
}

NumericVector dtmvnorm_botev(const arma::mat &x,
                             const arma::vec& mean,
                             const arma::mat& sigma,
                             const arma::vec& lower,
                             const arma::vec& upper){
  Environment env("package:TruncatedNormal"); 
  
  Function fun = env["dtmvnorm"];    
  
  // Convert arma to Rcpp quantities
  NumericMatrix x_rcpp(x.n_rows, x.n_cols, x.begin());
  NumericVector mu_rcpp(mean.begin(), mean.end());
  NumericMatrix sigma_rcpp(sigma.n_rows, sigma.n_cols, sigma.begin());
  NumericVector lb_rcpp(lower.begin(), lower.end());
  NumericVector ub_rcpp(upper.begin(), upper.end());
  
  NumericVector out = fun(Rcpp::_["x"]  = x_rcpp,
                          Rcpp::_["mu"] = mu_rcpp,
                          Rcpp::_["sigma"] = sigma_rcpp,
                          Rcpp::_["lb"] = lb_rcpp,
                          Rcpp::_["ub"] = ub_rcpp
                          ); 
  
  
  return out;
}

NumericVector dtmvnorm_botev(const arma::mat &x,
                             const arma::vec& mean,
                             const arma::mat& sigma,
                             const arma::vec& lower,
                             const arma::vec& upper,
                             bool indep){
  if(indep){
    Environment env("package:truncnorm"); 
    Function fun = env["dtruncnorm"];    
    
    NumericVector out(x.n_rows, 1.0);
    for(int j = 0; j < x.n_cols; ++j){
      NumericVector x_rcpp(x.col(j).begin(), x.col(j).end());
      NumericVector dens = fun(Rcpp::_["x"] = x_rcpp,
                               Rcpp::_["mean"] = mean(j),
                               Rcpp::_["sd"] = sqrt(sigma(j, j)),
                               Rcpp::_["a"] = lower(j),
                               Rcpp::_["b"] = upper(j));
      out = out * dens;
    }

    return out;
  } else {
    Environment env("package:TruncatedNormal"); 
    
    Function fun = env["dtmvnorm"];    
    
    // Convert arma to Rcpp quantities
    NumericMatrix x_rcpp(x.n_rows, x.n_cols, x.begin());
    NumericVector mu_rcpp(mean.begin(), mean.end());
    NumericMatrix sigma_rcpp(sigma.n_rows, sigma.n_cols, sigma.begin());
    NumericVector lb_rcpp(lower.begin(), lower.end());
    NumericVector ub_rcpp(upper.begin(), upper.end());
    
    NumericVector out = fun(Rcpp::_["x"]  = x_rcpp,
                            Rcpp::_["mu"] = mu_rcpp,
                            Rcpp::_["sigma"] = sigma_rcpp,
                            Rcpp::_["lb"] = lb_rcpp,
                            Rcpp::_["ub"] = ub_rcpp
    ); 
    return out;
  }
}
// //-------------------------------------------------------------------------
// //sample the labels 1:k with probability p 
// int sampleC(int k, arma::vec prob) {
// 	double rU;
// 	int j;
// 	double cumprob = 0;
// 	rU = unif_rand();
// 	double mass = 0;
// 	
// 	for (j = 0; j < k; j++) {
// 		mass = mass + prob(j);
// 	}
// 	for (j = 0; j < k; j++) {
// 		cumprob = cumprob + p(j)/mass;
// 		if (rU <= cumprob)
// 			break;
// 	}
// 	return j+1;
// }

//----------------------------------------------------------------------------
// get the scale probability to compare with the slice variable 
arma::vec scaleProb(const vector<double> &prob){
        int smax = n_scales_tree(prob);
        arma::vec pi(smax + 1, fill::zeros);
        for(int s = 0; s <= smax; s++) {
                int hmax = n_elem_scale(s); 
                for(int h = 1; h <= hmax; h++) 
                        pi(s) += extractNode(prob, s, h);
        }
        return pi;
}

//--------------------------------------------------------------------------
// Construct S and R trees [Equation 5 of Stefanucci and Canale (2021)]
void srTrees(
                vector<double> &Ssh,
                vector<double> &Rsh,
                const double &a,
                const double &b,
                const double &delta
        )
{
        if(Ssh.size() != Rsh.size()){
                throw "Ssh and Rsh trees are of unequal sizes";
        }
        int smax = n_scales_tree(Ssh);
        for(int s = 0; s <= smax; ++s){
                for(int h = 1; h <= n_elem_scale(s); ++h){
                        int idx = sh_to_idx(s, h);
                        if(s == smax){
                                Ssh[idx] = 1.0;
                                Rsh[idx] = 1.0;
                        } else {
                                Ssh[idx] = R::rbeta(1 - delta, a + delta*s);
                                Rsh[idx] = R::rbeta(b, b);
                        }
                }
        }
}


//-------------------------------------------------------------------------
// Compute probability tree from S and R trees
// [Equation 3 of Stefanucci and Canale (2021)]
vector<double> computeProb(
                vector<double> Stree,
                vector<double> Rtree,
                const double &a,
                const double &b)
{
    int smax = n_scales_tree(Stree);
    // Create auxiliary tree with root of Stree
    vector<double> tree(Stree.size(), 0.0);
    writeNode(tree, extractNode(Stree, 0, 1), 0, 1);
    
	int r, g_shr;
	int went_right;
	double T_shr;
	double I_S = 1;

	for(int s = 1; s <= smax; s++) {
		R_CheckUserInterrupt();
		for (int h = 1; h <= n_elem_scale(s); h++) {
			I_S = 1;
			for (r = 0; r < s; r++) {
				g_shr = ceil(h, (int) std::pow(2.0, s-r));
				went_right = ((2*g_shr) == ceil(h, (int)std::pow(2.0, (s-r-1))) ); //as above

				if (went_right) {
					T_shr = extractNode(Rtree, r, g_shr);
				}
				else {
					T_shr = 1 - extractNode(Rtree, r, g_shr);
				}

				I_S = I_S * (1 - extractNode(Stree, r, g_shr)) * T_shr;
			}
			writeNode(tree, extractNode(Stree, s, h)*I_S, s, h);
		}
	}
	return tree;
}

//-------------------------------------------------------------------------
// Construct auxiliary n,r,v trees given a cluster matrix
// All parameter trees (n, r, v) must be initialized to 0 in all nodes.
//
// @param n : empty binary tree containing the number of subjects assigned to node (s,h)
// @param r : binary tree containing the number of subjects which proceed to the right at node (s,h)
// @param v : binary tree containing the number of subjects which pass through node (s,h)
// @param cluster: Integer matrix of dimension n x 2 containing current cluster allocations
//
// @return none, modifies n, r, v inplace
void nrvTrees(
                vector<int> &n,
                vector<int> &r,
                vector<int> &v,
                const arma::umat &cluster
                )
{
        // Number of observations
        int N = cluster.n_rows;             
        int data_smax = cluster.col(0).max();

        // Cycle all nodes and construct n,r,v trees
        // Only access data up to maximum s in cluster, since everything below is zero
        for(int s = 0; s <= data_smax; ++s){
                checkUserInterrupt();
                int hmax = (int) std::pow(2.0, s);
                for(int h = 1; h <= hmax; ++h){
                        for(int i = 0; i < N; ++i){
                                int idx = sh_to_idx(s, h);
                                if((int) cluster(i, 0) == s && (int) cluster(i, 1) == h){
                                        n[idx] += 1;
                                        v[idx] += 1;
                                }
                                int rel_s = cluster(i, 0) - s;
                                if (rel_s > 0) {
                                        if( ceil( cluster(i, 1), ( (int) pow(2.0, rel_s))) == h) v[idx] += 1;
                                        if( ceil( cluster(i, 1), ( (int) pow(2.0, rel_s - 1))) == 2 * h) r[idx] += 1;
                                }
                        }
                }
        }
}

//' Compute prior distribution of Multiscale Mixture of Kernels
//'
//' @param TH: Binary tree of location parameters
//' @param SIG: Binary tree of scale matrices
//' @param prob: Binary tree of node probabilities
//' @param Ssh: Binary tree of conditional stopping probabilities
//' @param Rsh: Binary tree of conditional right-turn probabilities
//' @param thrs: Binary tree of matrices with thresholds for each dimension
//' @param mu0: Prior mean vector of G_0.
//' @param k0: Prior marginal variances of G_0.
//' @param sig0: Prior covariance matrix of H_0.
//' @param a: Prior stick-breaking weight.
//' @param b: Prior stick-breaking weight.
//' @param delta: Prior stick-breaking weight.
//' @param indep: If TRUE assumes diagonal scale matrix.
//' @param verbose: If TRUE prints information about the constructed variables
//'
//' @return List of values pertaining to the prior multiscale distribution
void msMK_prior(arma::mat &TH,
                arma::cube &SIG,
                const arma::cube &thrs,
                const arma::colvec& mu0,
                const arma::colvec& k0,
                const arma::mat& sig0,
                const double& a,
                const double& b,
                const double& delta,
                const bool& indep
        )
{
        int p = sig0.n_cols;

        int smax = n_scales_tree(TH);
        for(int s = 0; s <= smax; ++s){
               for(int h = 1; h <= n_elem_scale(s); ++h){
                 mat thrs_sh = extractNode(thrs, s, h); 
                 // Generate locations from prior location process
                 vec theta(p);
                 for(int j = 0; j < p; ++j){
                   theta(j) = rtmvnorm_arma_botev(mu0(j), sqrt(k0(j)), thrs_sh(0, j), thrs_sh(1, j));
                 }
                 writeNode(TH, theta, s, h);
                 
                 // Generate scale matrices from prior scale process
                 mat sigma(p, p);
                 if(indep){
                   // Inverse gamma (diagonal matrices)
                   vec scale(p);
                   for(int j = 0; j < p; ++j){
                     scale(j) = rinvgamma(std::pow(2.0, 6), std::pow(2.0, 6 - s) * sig0(j, j));
                   }
                   sigma = arma::diagmat(scale);
                 } else {
                   // Inverse wishart
                   sigma = riwish((int) pow(2.0, s) + p + 1, sig0);
                 }
                 writeNode(SIG, sigma, s, h);
               } 
        }
}

arma::vec rowSum(const arma::mat &x){
    int n = x.n_rows;
    int p = x.n_cols;
    vec out(n, fill::zeros);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < p; ++j){
            out(i) += x(i,j);
        }
    }
    return out;
}
