#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of similarity matrix computation using precomputed distances
//'
//'
//'@param c an MCMC partitions of length \code{n}.
//'
//'@param d a symmetric \code{n x n} matrix containing distances
//'between each group distributions.
//'
//'@author Boris Hejblum, Chariff Alkhassim
//'
//'@export
//'
//'@examples
//'c <- c(1,1,2,3,2,3)
//'d <- matrix(runif(length(c)^2),length(c))
//'NuMatParC(c,d)
//'
//'
// [[Rcpp::export]]
List NuMatParC(NumericVector c, arma::mat d){

  int n = c.size();

  mat similarity = mat(n, n, fill::zeros);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      similarity(i,j)=(c(i) != c(j) )*d(c(i)-1,c(j)-1);
      similarity(j,i) = similarity(i,j);
    }
  }
return Rcpp::List::create(Rcpp::Named("NuMatParC") = similarity);
}
