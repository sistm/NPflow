#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Parallel C++ implementation of the multinomial sampling from a matrix 
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
IntegerVector sampleClassC_bis_par(NumericMatrix probMat,
                                   int ncores=1){
    
    const mat pm =as<mat>(probMat);
    const int n =pm.n_cols;
    const int p = pm.n_rows;
    IntegerVector c = IntegerVector(n);
        
    omp_set_num_threads(ncores);
    
    #pragma omp parallel for shared(c)    
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
