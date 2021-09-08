/*--------------------------------------------------------------------------
 Multiscale Mixtures of Kernels [msMK]
 gibbs.cpp - C++ functions to perform Gibbs sampling for the multiscale mixture of Gaussian kernels
 
 Version 0.1 of September 2021
 2013 - Antonio Canale (antonio.canale@unipd.it)
 2021 - Daniele Zago (daniele.zago.1@studenti.unipd.it)
 --------------------------------------------------------------------------*/
#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <vector>
#include "bintree.h"
#include "auxGibbs.h"
#include "gibbs.h"

using namespace arma;
using namespace Rcpp;
using namespace std;
using std::vector;

//' Allocate each observation to a tree cluster, conditionally on the current
//' value oflocation, scale, and probabilities -- [Algorithm 2]
//'
//' @param TH: field of location parameters
//' @param SIG: field of scale matrices 
//' @param thrs: field of matrices containing lower and upper truncation bounds
//' @param prob: vector of node probability weights
//' @param pi_s: vector of scale probabilities such that pi_s(i) = sum_h pi_sh(i, h)
//' @param data: matrix observed data
//' @param indep: bool indicating whether kernels are independent or full-covariance normals
//'
//' @return two-column matrix where every row (s_i, h_i) is the cluster allocation of the i-th data point
arma::umat clusterAllocation(
                const arma::mat &TH,
                const arma::cube &SIG,
                const arma::cube &thrs,
                const vector<double> &prob,
                const arma::vec &pi_s,
                const arma::mat &data,
                bool indep
        )
{
        const int n = data.n_rows;

        // Matrix that will be returned by the function
        umat cl_new(n, 2, fill::zeros);

        int smax = n_scales_tree(TH);
        int nelem = smax_to_nelem(smax);
        
        mat dens_mat(n, nelem, fill::zeros);
        mat post_s(n, smax + 1, fill::zeros);
        double pi_old = pi_s.max();

        for(int s = 0; s <= smax; ++s){
            int hmax = n_elem_scale(s);
            for(int h = 1; h <= hmax; ++h){
                vec mu_sh = extractNode(TH, s, h);
                // Matrix of thresholds
                mat sigma_sh = extractNode(SIG, s, h);
                mat thrs_sh = extractNode(thrs, s, h);  
                // Lower truncation bounds
                vec lb = thrs_sh.row(0).t();            
                // Upper truncation bound
                vec ub = thrs_sh.row(1).t();            
                // Tree to vector index
                int idx = sh_to_idx(s, h);              
                // Calculate conditional weight for scale s
                double barpi_sh = prob[idx] / pi_s(s); 

                // Density at scale (s,h)
                vec dens = dmvnorm(data, mu_sh, sigma_sh);
                dens_mat.col(idx) = dens;
                post_s.col(s) += barpi_sh * dens;
            }
            
            // post_s(s) is zero if u > pi_s(s)
            for(int i = 0; i < n; ++i){
                double u = R::runif(0.0, pi_old);
                if(u > pi_s(s)){
                    for(int h = 1; h <= hmax; ++h){
                        post_s(i, s) = 0.0;
                    }
                }
            }
        }
        
        // Normalize probabilities
        vec rsums = rowSum(post_s);
        for(int i = 0; i < n; ++i){
            rowvec post_s_i = post_s.row(i) / rsums(i);
            int s_i = Rcpp::RcppArmadillo::sample(linspace<uvec>(0, smax, smax + 1), 1, FALSE, post_s_i.t())(0);
            int hs = n_elem_scale(s_i);
            int idx_i_start = (int) pow(2.0, s_i) - 1;
            int idx_i_end = (int) std::pow(2.0, s_i) - 2 + hs;
            
            // Extract probabilities for given scale
            vec arma_prob = as<arma::colvec>(wrap(prob));
            
            vec barpi_i = arma_prob.subvec(idx_i_start, idx_i_end) / pi_s(s_i);
            
            // Extract densities for given scale
            vec dens_i = (dens_mat.row(i)).t();
            vec dens_i_sh = dens_i.subvec(idx_i_start, idx_i_end);
            
            // Vector of probabilities for node h conditional on scale s
            // Equation (12), elementwise product
            vec post_h = barpi_i % dens_i_sh;
            post_h /= accu(post_h);
            
            // Sample scales {1, 2, ..., hs} according to probability post_h
            int h_i = Rcpp::RcppArmadillo::sample(linspace<uvec>(1, hs, hs), 1, FALSE, post_h)(0);
            cl_new(i, 0) = s_i;
            cl_new(i, 1) = h_i;
        }
        
        return cl_new;
}

//' Allocate each observation to a tree cluster, conditionally on the current
//' value of location, scale, and probabilities given previous cluster allocation
//' -- [Algorithm 2]
//'
//' @param TH: field of location parameters
//' @param SIG: field of scale matrices 
//' @param thrs: field of matrices containing lower and upper truncation bounds
//' @param prob: vector of node probability weights
//' @param pi_s: vector of scale probabilities such that pi_s(i) = sum_h pi_sh(i, h)
//' @param data: matrix observed data
//' @param cl_old: previous cluster assignment matrix
//' @param indep: bool indicating whether kernels are independent or full-covariance normals
//'
//' @return two-column matrix where every row (s_i, h_i) is the cluster allocation of the i-th data point
arma::umat clusterAllocation(
                const arma::mat &TH,
                const arma::cube &SIG,
                const arma::cube &thrs,
                const vector<double> &prob,
                const arma::vec &pi_s,
                const arma::mat &data,
                const arma::umat &cl_old,
                bool indep
        )
{
        const int n = data.n_rows;

        // Matrix that will be returned by the function
        umat cl_new(n, 2, fill::zeros);

        int smax = n_scales_tree(TH);
        int nelem = smax_to_nelem(smax);
        
        mat dens_mat(n, nelem, fill::zeros);
        mat post_s(n, smax + 1, fill::zeros);
        
        for(int s = 0; s <= smax; ++s){
            int hmax = n_elem_scale(s);
            for(int h = 1; h <= hmax; ++h){
                vec mu_sh = extractNode(TH, s, h);
                mat sigma_sh = extractNode(SIG, s, h);
                // Matrix of thresholds
                mat thrs_sh = extractNode(thrs, s, h);  
                // Lower truncation bounds
                vec lb = thrs_sh.row(0).t();            
                // Upper truncation bound
                vec ub = thrs_sh.row(1).t();            
                // Tree to vector index
                int idx = sh_to_idx(s, h);              

                // Calculate conditional weight for scale s
                double barpi_sh = prob[idx] / pi_s(s); 

                // Density at scale (s,h)
                vec dens = dmvnorm(data, mu_sh, sigma_sh);
                dens_mat.col(idx)  = dens;
                post_s.col(s) += barpi_sh * dens;
            }
            
            // post_s(s) is zero if u > pi_s(s)
            for(int i = 0; i < n; ++i){
                double u = R::runif(0.0, pi_s(cl_old(i, 0)));
                if(u > pi_s(s)){
                    for(int h = 1; h <= hmax; ++h){
                        post_s(i, s) = 0.0;
                    }
                }
            }
        }
        
        // Normalize probabilities
        vec rsums = rowSum(post_s);
        for(int i = 0; i < n; ++i){
            rowvec post_s_i = post_s.row(i) / rsums(i);
            int s_i = Rcpp::RcppArmadillo::sample(linspace<uvec>(0, smax, smax + 1), 1, FALSE, post_s_i.t())(0);
            int hs = n_elem_scale(s_i);
            int idx_i_start = (int) pow(2.0, s_i) - 1;
            int idx_i_end = (int) std::pow(2.0, s_i) - 2 + hs;
            
            // Extract probabilities for given scale
            vec arma_prob = as<arma::colvec>(wrap(prob));
            vec barpi_i = arma_prob.subvec(idx_i_start, idx_i_end) / pi_s(s_i);
            
            // Extract densities for given scale
            vec dens_i = (dens_mat.row(i)).t();
            vec dens_i_sh = dens_i.subvec(idx_i_start, idx_i_end);
            
            // Vector of probabilities for node h conditional on scale s
            // Equation (12), elementwise product
            vec post_h = barpi_i % dens_i_sh;
            post_h /= accu(post_h);
            
            // Sample scales {1, 2, ..., hs} according to probability post_h
            int h_i = Rcpp::RcppArmadillo::sample(linspace<uvec>(1, hs, hs), 1, FALSE, post_h)(0);
            cl_new(i, 0) = s_i;
            cl_new(i, 1) = h_i;
        }
        
        return cl_new;
}

//' Update tree probability weights conditional on cluster allocation 
//' (Algorithm 3 of Canale and Dunson, 2014)
//'
//' @param cluster: matrix of dimension n x 2 containing current cluster allocations
//' @param smax: int specifying maximum tree depth 
//' @param a: stick-breaking parameter
//' @param b: stick-breaking parameter
//' @param delta: stick-breaking parameter
//'
//' @return tree containing updated probabilities for each node
std::vector<double> treeUpdate(
        const arma::umat& cluster,
        const int &smax,
        const double &a,
        const double &b,
        const double &delta
        )
{
    // Initialize n,r,v trees based on current cluster allocations
    int nelem = smax_to_nelem(smax);
    vector<int> n(nelem, 0);
    vector<int> r(nelem, 0);
    vector<int> v(nelem, 0);
    nrvTrees(n, r, v, cluster);
    
    // Initialize S and R trees
    vector<double> Stree(nelem, 0.0);
    vector<double> Rtree(nelem, 0.0);
    
    for(int s = 0; s <= smax; ++s){
        // Number of elements in scale
        int hmax = n_elem_scale(s);
        
        // Update Ssh
        if(s != smax){
            // Generate from Beta(1 - delta*n_sh, a + delta(s+1) + v_sh - n_sh)
            for(int h = 1; h <= hmax; ++h){
                double Stree_sh = R::rbeta(1.0 - delta + (double) extractNode(n, s, h),
                                       a + delta * ( (double) s + 1.0) +
                                           (double) extractNode(v, s, h) -
                                           (double) extractNode(n, s, h));
                writeNode(Stree, Stree_sh, s ,h);
            }
        } else{
            // Pad S with ones if it is the last level
            for(int h = 1; h <= hmax; ++h){
                writeNode(Stree, 1.0, s, h);
            }
        }
        // Update Rsh from Beta(b + r_sh, b + v_sh - n_sh - r_sh)
        for(int h = 1; h <= hmax; ++h){
            double Rtree_sh = R::rbeta(b + (double) extractNode(r, s, h),
                 b + (double) extractNode(v, s, h) - (double) extractNode(n, s, h)
                                   - (double) extractNode(r, s, h));
            writeNode(Rtree, Rtree_sh, s, h);
        }
    }

   // Compute tree of probabilities
   return computeProb(Stree, Rtree, a, b);
}




//' Gibbs sampler update of location parameters based on current cluster allocation and value of scale parameters
//'
//' @param SIG: field of scale matrices 
//' @param thrs: field of matrices containing lower and upper truncation bounds
//' @param prob: vector of node probability weights
//' @param mu0: prior location parameter
//' @param sig0: prior location matrix
//' @param cluster: matrix of integer of size n_obs x 2 containing current cluster allocations
//' @param data: matrix of observed data 
//' @param indep: bool indicating whether kernels are independent or full-covariance normals
//'
//' @return List containing the updated location parameters
arma::mat thetaUpdate(
        const arma::cube &SIG,
        const arma::cube &thrs,
        const arma::mat &mu0,
        const arma::mat &sig0,
        const arma::umat &cluster,
        const arma::mat &data,
        bool indep
    )
{
    int p = data.n_cols;
    // Maximum depth
    int smax = n_scales_tree(SIG);   
    // Number of elements in array representation
    int nelem = smax_to_nelem(smax); 
    // Inverse of sig0
    mat sig0_i = sig0.i();           
    
    mat thetaNew(p, nelem);

    for(int s = 0; s <= smax; ++s){
        int hmax = n_elem_scale(s);

        for(int h = 1; h <= hmax; ++h){
            // Container for sampled location
            vec theta(p, fill::zeros);              
            // Inverse of scale at node (s,h)
            mat sigma = extractNode(SIG, s, h);     
            // Matrix of thresholds 
            mat thrs_sh = extractNode(thrs, s, h);  
            // Lower truncation bounds
            vec lb = thrs_sh.row(0).t();           
            // Upper truncation bound
            vec ub = thrs_sh.row(1).t();           

            // Get currently allocated observations
            uvec data_idx = find((cluster.col(0) == s) && (cluster.col(1) == h));
            mat data_now = data.rows(data_idx);
            int n_now = data_now.n_rows;

            if(n_now != 0){
            // If observations are allocated, update using Gibbs sampler formulas
                if(indep){
                    vec ysum = n_now * mean(data_now).t();
                    for(int j = 0; j < p; ++j){
                        double new_s = sigma(j, j)*sig0(j, j) / (n_now * sig0(j, j) + sigma(j, j));
                        double new_m = (mu0(j) * sigma(j,j) + ysum(j) * sig0(j,j))/ (n_now * sig0(j,j) + sigma(j,j));
                        theta(j) = rtmvnorm_arma_botev(new_m, new_s, lb(j), ub(j));
                    }
                } else {
                    // Column mean
                    vec ybar = mean(data_now).t(); 
                    mat sigma_i = inv_sympd(sigma);         
                    mat S_n = (sig0_i + n_now * sigma_i).i();
                    vec mu_n = S_n * (sig0_i * mu0 + n_now * sigma_i * ybar);
                    theta = rtmvnorm_arma_botev(mu_n, S_n, lb, ub);
                }
            } else {
                // If no observations are allocated, use prior distribution
                if(indep){
                   for(int j = 0; j < p; ++j){
                        theta(j) = rtmvnorm_arma_botev(mu0(j), sig0(j, j), lb(j), ub(j));
                   }

                } else {
                    theta = rtmvnorm_arma_botev(mu0, sig0, lb, ub);
                }
            }
            writeNode(thetaNew, theta, s, h);
        }
        
    }
    return thetaNew;
}


//' Gibbs sampler update of scale parameters based on current cluster allocation and value of location parameters
//'
//' @param TH: field of location parameters 
//' @param mu0: prior location parameter
//' @param sig0: prior location matrix
//' @param cluster: matrix of integer of size n_obs x 2 containing current cluster allocations
//' @param data: matrix of observed data 
//' @param indep: bool indicating whether kernels are independent or full-covariance normals
//'
//' @return List containing the updated scale parameters
arma::cube sigmaUpdate(
        const arma::mat TH,
        const arma::mat &mu0,
        const arma::mat &sig0,
        const arma::umat &cluster,
        const arma::mat &data,
        bool indep)
{
    
    int p = data.n_cols;
    // Maximum depth
    int smax = n_scales_tree(TH);   
    // Number of elements in array representation
    int nelem = smax_to_nelem(smax); 

    cube sigmaNew(p, p, nelem);

    for(int s = 0; s <= smax; ++s){
        int hmax = n_elem_scale(s);

        for(int h = 1; h <= hmax; ++h){
            // Container for sampled scale matrix
            mat sigma(p, p, fill::zeros);           
            // Current location parameter
            vec mu_sh = extractNode(TH, s, h);      
            
            // Get currently allocated observations
            uvec idx = find((cluster.col(0) == s) && (cluster.col(1) == h));
            mat data_now = data.rows(idx);
            int n_now = data_now.n_rows;

            if(n_now != 0){
                // If observations are allocated, update using Gibbs sampler formulas
                if(indep){
                    // Independent updates
                    double par1 = std::pow(2.0, 6) + n_now/2;
                    double par2;
                    for(int j = 0; j < p; ++j){
                        vec diff = data_now.col(j) - mu_sh(j);
                        par2 = std::pow(2.0, 6 - s) * sig0(j,j) + accu(diff % diff) / 2;
                        sigma(j, j) = rinvgamma(par1, par2);
                    }
                } else {
                    // Dependent updates
                    int nu_n = 1 + p + std::pow(2.0, s) + n_now;
                    // Matrix of sum of squares
                    mat sig_n(p, p, fill::zeros);                       
                    for(int i = 0; i < n_now; ++i){
                        vec v = ((data_now.row(i)).t() - mu_sh);
                        sig_n += v * v.t();
                    }
                    mat psi_n = sig0 + sig_n;
                    sigma = riwish(nu_n, psi_n);
                }
            } else {
                // If no observations are allocated, use prior distribution
                if(indep){
                    // Independent updates
                    for(int j = 0; j < p; ++j){
                        sigma(j, j) = rinvgamma(std::pow(2.0, 6), std::pow(2.0, 6 - s) * sig0(j, j));
                    }
                } else {
                    // Dependent updates
                    sigma = riwish(std::pow(2.0, s) + p + 1, sig0);
                }
            }
            writeNode(sigmaNew, sigma, s, h);
        }
    }
    return sigmaNew;
}
