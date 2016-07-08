#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of the multinomial sampling from a matrix
//' of column vectors, each containing the sampling probabilities
//' for their respective draw
//'
//'@param probMat a numeric matrix of dim \code{k x n} of containing column vectors of sampling
//'probabilities for each class \code{k}.
//'
//'@return a vector of integer of length \code{n} containing the multinomial draws for each
//'observation, i.e. the class allocation.
//'
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
IntegerVector sampleClassC(arma::mat probMat){

    int n = probMat.n_cols;
    IntegerVector c(n);

    GetRNGstate();

    for(int i=0; i<n; i++){
        vec pm_i = probMat.col(i);
        vec probs = pm_i/sum(pm_i);
        vec probs_cumsum = cumsum(probs);
        uvec c_i = find(probs_cumsum > unif_rand());
        c(i) = c_i(0) + 1;
    }

    PutRNGstate();

    return(c);
}
