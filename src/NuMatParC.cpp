#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation
//' 
//'
//'@param \code{c} is an MCMC partitions of length \code{n}.
//'
//'@param \code{d} is a symmetric \code{n x n} matrix containing distances
//'between each group distributions.
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
List NuMatParC(NumericVector c, NumericMatrix d){
  
  int n = c.size();
  mat dd = as<mat>(d);
  
  mat similarity = mat(n, n, fill::zeros);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      similarity(i,j)=(c(i) != c(j) )*dd(c(i)-1,c(j)-1);
      similarity(j,i) = similarity(i,j);
    }
  }
return Rcpp::List::create(Rcpp::Named("NuMatParC") = similarity);
}
