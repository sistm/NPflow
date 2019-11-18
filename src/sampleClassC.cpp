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
//'@param Log a logical flag indicating whether the provided \code{probMat} is on the log scale
//'or natural probability scale. Default is \code{FALSE} in which case it is considered on the natural
//'probability scale.
//'
//'@return a vector of integer of length \code{n} containing the multinomial draws for each
//'observation, i.e. the class allocation.
//'
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
IntegerVector sampleClassC(const arma::mat & probMat,
                           const bool & Log=false){

    int n = probMat.n_cols;
    IntegerVector c(n);

    GetRNGstate(); //initialize Random Number Generator state
    if(Log){
      for(int i=0; i<n; i++){
        vec pm_i = probMat.col(i);
        double pmax_i = max(pm_i);
        vec probs =  pm_i - (pmax_i + log(sum(exp(pm_i - pmax_i))));
        double probsmax = max(probs);
        vec probs_cumsum = exp(probsmax + log(cumsum(exp(probs-probsmax))));
        uvec c_i = find(probs_cumsum > unif_rand());
        c(i) = c_i(0) + 1;
      }
    }else{
      for(int i=0; i<n; i++){
          vec pm_i = probMat.col(i);
          vec probs = pm_i/sum(pm_i);
          vec probs_cumsum = cumsum(probs);
          uvec c_i = find(probs_cumsum > unif_rand());
          c(i) = c_i(0) + 1;
      }
    }
    PutRNGstate(); //export Random Number Generator state

    return(c);
}
