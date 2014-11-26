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
//'@param probMat
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
IntegerVector sampleClassC_par(NumericMatrix probMat,
                               int ncores=1){
    
    const mat pm =as<mat>(probMat);
    const int n = pm.n_cols;
    IntegerVector c = IntegerVector(n);
    
    omp_set_num_threads(ncores);
    
    #pragma omp parallel for shared(c)    
    for(int i=0; i<n; i++){
        vec pm_i = pm.col(i);
        vec probs = pm_i/sum(pm_i);
        vec probs_cumsum = cumsum(probs);
        uvec c_i = find(probs_cumsum > unif_rand());
        c(i) = c_i(0) + 1;  
    }
    
    return(c);
}
