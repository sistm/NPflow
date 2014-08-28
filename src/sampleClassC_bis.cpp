#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of the multinomial sampling from a matrix 
//' of column vectors each containing the sampling probabilities 
//' for their respective draw
//' 
//' @details Slower than sampleClassC
//' 
//'@param probMat
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
IntegerVector sampleClassC_bis(NumericMatrix probMat){
    
    mat pm =as<mat>(probMat);
    int n =pm.n_cols;
    int p = pm.n_rows;
    IntegerVector c = IntegerVector(n);
    
    for(int i=0; i<n; i++){
        vec pm_i = pm.col(i);
        NumericVector probs = as<NumericVector>(wrap(pm_i));
        probs = probs/sum(probs);
        IntegerVector ans = IntegerVector(p);
        R::rmultinom(1, probs.begin(), p, ans.begin());
        NumericVector s = as<NumericVector>(wrap(ans));
        vec s2 = as<vec>(s);
        uvec c_i = find(s2);
        c(i) = c_i(0) + 1;
        
    }
    
    return(c);
}
