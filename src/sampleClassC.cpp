#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of the multinomial sampling from a matrix 
//' of column vectors each containing the sampling probabilities 
//' for their respective draw
//' 
//'@param probMat
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
IntegerVector sampleClassC(NumericMatrix probMat){
    
    mat pm =as<mat>(probMat);
    int n = pm.n_cols;
    IntegerVector c = IntegerVector(n);
    
    for(int i=0; i<n; i++){
        vec pm_i = pm.col(i);
        vec probs = pm_i/sum(pm_i);
        vec probs_cumsum = cumsum(probs);
        uvec c_i = find(probs_cumsum > unif_rand());
        c(i) = c_i(0) + 1;
        
    }
    
    return(c);
}
