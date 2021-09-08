/*--------------------------------------------------------------------------
 Multiscale Mixtures of Kernels [msMK]
 Bintree.h - Auxiliary C++ functions to deal with binary tree structures
 
 Version 0.1 of September 2021 
 2013 - Antonio Canale (antonio.canale@unipd.it)
 2021 - Daniele Zago (daniele.zago.1@studenti.unipd.it)
 --------------------------------------------------------------------------*/
#include <math.h>
#include <R.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <math.h>

using namespace std;

#ifndef BINTREE_H_
#define BINTREE_H_

inline int sh_to_idx(int s, int h){
    return (int) std::pow(2.0, s) - 2 + h;
}

inline int smax_to_nelem(int smax){
    return (int) std::pow(2.0, smax + 1) - 1;
}

inline int n_elem_scale(int s){
    return (int) std::pow(2.0, s);
}


// Return number of scales (s) in the tree ----------------------
template<typename T>
int n_scales_tree(const vector<T> &tree){
    return ((int) std::log2(tree.size() + 1)) - 1;
        
}

inline int n_scales_tree(const arma::vec &tree){
    return ((int) std::log2(tree.n_elem + 1)) - 1;
        
}

inline int n_scales_tree(const arma::mat &tree){
    return ((int) std::log2(tree.n_cols + 1)) - 1;
        
}

inline int n_scales_tree(const arma::cube &tree){
    return ((int) std::log2(tree.n_slices + 1)) - 1;
        
}


// Write node to a location in the tree ----------------------
template<typename T>
void writeNode(vector<T> &tree, const T &node,  const int &s, const int &h){
    tree[sh_to_idx(s, h)] = node;
}

inline void writeNode(arma::vec &tree, const double &node,  const int &s, const int &h){
    tree(sh_to_idx(s, h)) = node;
}

inline void writeNode(arma::mat &tree, const arma::colvec &node,  const int &s, const int &h){
    tree.col(sh_to_idx(s, h)) = node;
}

inline void writeNode(arma::cube &tree, const arma::mat &node,  const int &s, const int &h){
    tree.slice(sh_to_idx(s, h)) = node;
}

// Extract node from a location in the tree ----------------------
template<typename T>
T extractNode(const vector<T> &tree, const int &s, const int &h){
    return tree[sh_to_idx(s, h)];
}

inline double extractNode(const arma::vec &tree, const int &s, const int &h){
  return tree(sh_to_idx(s, h));
}

inline arma::colvec extractNode(const arma::mat &tree, const int &s, const int &h){
    return tree.col(sh_to_idx(s, h));
}

inline arma::mat extractNode(const arma::cube &tree, const int &s, const int &h){
    return tree.slice(sh_to_idx(s, h));
}


// Print contents of the tree to standard output ----------------------
template<typename T>
void printTree(const vector<T> &tree){
    int smax = n_scales_tree(tree);
    
    for(int s = 0; s <= smax; ++s){
        Rcpp::Rcout << "-- s = " << s << " --" << endl;
        for(int h = 1; h <= (int) std::pow(2.0, s); ++h){
            Rcpp::Rcout << extractNode(tree, s, h) << endl; 
        }
    }
    Rcpp::Rcout << endl;
}

inline void printTree(const arma::vec &tree){
    int smax = n_scales_tree(tree);
    
    for(int s = 0; s <= smax; ++s){
        Rcpp::Rcout << "-- s = " << s << " --" << endl;
        for(int h = 1; h <= (int) std::pow(2.0, s); ++h){
            Rcpp::Rcout << extractNode(tree, s, h) << endl; 
        }
    }
    Rcpp::Rcout << endl;
}

inline void printTree(const arma::mat &tree){
    int smax = n_scales_tree(tree);
    
    for(int s = 0; s <= smax; ++s){
        Rcpp::Rcout << "-- s = " << s << " --" << endl;
        for(int h = 1; h <= (int) std::pow(2.0, s); ++h){
            Rcpp::Rcout << extractNode(tree, s, h) << endl; 
        }
    }
    Rcpp::Rcout << endl;
}

inline void printTree(const arma::cube &tree){
    int smax = n_scales_tree(tree);
    
    for(int s = 0; s <= smax; ++s){
        Rcpp::Rcout << "-- s = " << s << " --" << endl;
        for(int h = 1; h <= (int) std::pow(2.0, s); ++h){
            Rcpp::Rcout << extractNode(tree, s, h) << endl; 
        }
    }
    Rcpp::Rcout << endl;
}


// Conditions for the 
inline int compute_right_h(int s, int h) {
  int right_h = h - (int) std::pow(2.0, s - 1);
  return right_h;
}

inline double compute_condition(int s, int h) {
  double h_ = (double)h;
  double s_ = (double)s;
  double cond = h_ / std::pow(2.0, s_);
  return cond;
}

#endif /* BINTREE_H_ */
