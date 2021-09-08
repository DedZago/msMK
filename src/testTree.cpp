#include <RcppArmadillo.h>
#include "bintree.h"
#include "auxGibbs.h"
#include "gibbs.h"
#include <vector>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]

// [[Rcpp::export]]
Rcpp::List msMK_mcmc_test(
    const int& sample,
    const arma::mat &y,
    const int &a,
    const int &b,
    const double &delta,
    const int &smax,
    const arma::vec &mu0,
    const arma::vec &k0,
    const arma::mat &sig0,
    bool indep,
    const arma::mat &lbNumpy,
    const arma::mat &ubNumpy,
    const int &burnin
)
{
  int nsim = sample + burnin;
  int p = mu0.n_elem;
  int nelem = smax_to_nelem(smax);
  
  // Initialize Stree, Rtree
  vector<double> Stree(nelem, 0.0);
  vector<double> Rtree(nelem, 0.0);
  
  // Rcout << "srTrees()\n" << endl;
  
  // Fill srTrees
  srTrees(Stree, Rtree, a, b, delta);
  
  // Calculate probability tree
  vector<double> prob = computeProb(Stree, Rtree, a, b);
  
  // Store thresholds from numpy matrices
  cube thrs(p, p, nelem);
  
  for(int i = 0; i < nelem; ++i){
    mat temp_thrs(2, p, fill::zeros);
    temp_thrs.row(0) = lbNumpy.row(i); 
    // Rcout << lbNumpy.row(i);
    temp_thrs.row(1) = ubNumpy.row(i); 
    // Rcout << ubNumpy.row(i);
    thrs.slice(i) = temp_thrs;
  }

  mat TH(p, nelem);
  cube SIG(p, p, nelem);

  // Rcout << "msMK_prior()\n" << endl;
  msMK_prior(TH, SIG, thrs, mu0, k0, sig0, a, b, delta, indep);
  // printTree(TH);
  // printTree(SIG);

  // Rcout << "clusterAllocation()\n" << endl;
  vec pi_s = scaleProb(prob);
  umat cluster = clusterAllocation(TH, SIG, thrs, prob, pi_s, y, indep);

  // Rcout << "treeUpdate()\n" << endl;
  prob = treeUpdate(cluster, smax, a, b, delta);
  // Rcout << "prob\n";
  // printTree(prob);

  // Rcout << "thetaUpdate()\n" << endl;
  mat thetaNew = thetaUpdate(SIG, thrs, mu0, sig0, cluster, y, indep);
  // Rcout << "thetaNew\n";
  // printTree(thetaNew);
  
  // Rcout << "sigmaUpdate()\n" << endl;
  cube sigmaNew = sigmaUpdate(TH, mu0, sig0, cluster, y, indep);

  // Lists of MCMC outputs
  field<mat> outTH(nsim);
  field<cube> outSIG(nsim);
  ucube outCluster(y.n_rows, 2, nsim);
  field<vector<double> > outProb(nsim);
  
  Rcout << "Starting iterations" << endl;
  // Apply MCMC
  for(int i = 0; i < nsim; ++i){
    if((i % 100) == 0){
      if(i < burnin){
        Rcout << "Burnin -- " << i << "/" << burnin << endl;
        if(i == burnin){
          Rcout << "Sampling -- " << 1 << "/" << nsim - burnin << endl;
        }
      }
      else{
        Rcout << "Sampling -- " << i - burnin << "/" << nsim - burnin << endl;
      }
    }
    pi_s = scaleProb(prob);
    // Rcout << "pi_s:\n";
    // Rcout << pi_s << endl;
    cluster = clusterAllocation(thetaNew, sigmaNew, thrs, prob, pi_s, y, cluster, indep);
    // Rcout << "cluster:\n";
    // Rcout << cluster << endl;
    outCluster.slice(i) = cluster;
    prob = treeUpdate(cluster, smax, a, b, delta);
    // Rcout << "prob:\n";
    // printTree(prob);
    outProb(i) = prob;
    thetaNew = thetaUpdate(sigmaNew, thrs, mu0, sig0, cluster, y, indep);
    // Rcout << "thetaNew:\n";
    // printTree(thetaNew);
    outTH(i) = thetaNew;
    sigmaNew = sigmaUpdate(thetaNew, mu0, sig0, cluster, y, indep);
    // Rcout << "sigmaNew:\n";
    // printTree(sigmaNew);
    outSIG(i) = sigmaNew;
  }
  Rcout << "Sampling -- " << nsim - burnin << "/" << nsim - burnin << endl;
  
  return Rcpp::List::create(
    Named("theta") = outTH,
    Named("sigma") = outSIG,
    Named("prob")  = outProb,
    Named("cluster") = outCluster,
    Named("thrs")  = thrs,
    Named("burnin") = burnin,
    Named("indep") = indep
    );
}

