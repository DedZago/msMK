#include <math.h>
#include <R.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <vector>
#include "bintree.h"
#include "auxGibbs.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef GIBBS_H_ 
#define GIBBS_H_

arma::umat clusterAllocation(
                const arma::mat &TH,
                const arma::cube &SIG,
                const vector<double> &prob,
                const arma::vec &pi_s,
                const arma::mat &data,
                bool indep);

arma::umat clusterAllocation(
                const arma::mat &TH,
                const arma::cube &SIG,
                const vector<double> &prob,
                const arma::vec &pi_s,
                const arma::mat &data,
                const arma::umat &cl_old,
                bool indep,
                NumericVector &cpo_inv
                );

vector<double> treeUpdate(
        const arma::umat& cluster,
        const int &smax,
        const double &a,
        const double &b,
        const double &delta);

arma::mat thetaUpdate(
        const arma::cube &SIG,
        const arma::cube &thrs,
        const arma::mat &mu0,
        const arma::mat &sig0,
        const arma::umat &cluster,
        const arma::mat &data,
        bool indep);

arma::cube sigmaUpdate(
        const arma::mat TH,
        const arma::mat &mu0,
        const arma::mat &sig0,
        const arma::umat &cluster,
        const arma::mat &data,
        bool indep);
        
# endif /* GIBBS_H_ */