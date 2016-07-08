#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate normal likelihood function for multiple inputs
//'
//'Based on the implementation from Nino Hardt and Dicko Ahmadou
//'http://gallery.rcpp.org/articles/dmvnorm_arma/
//'(accessed in August 2014)
//'
//'@param x data matrix
//'@param mean mean vector
//'@param varcovM variance covariance matrix
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}
//'@return vector of densities
//'
//'@author Boris P. Hejblum
//'
//'@export
//'
//'@examples
//'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
//'mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
//'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
//'mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1))
//'
//'if(require(microbenchmark)){
//' library(microbenchmark)
//' microbenchmark(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                times=10000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
//'
// [[Rcpp::export]]
NumericVector mvnpdfC(NumericMatrix x,
                      NumericVector mean,
                      NumericMatrix varcovM,
                      bool Log=true){

  mat xx = as<mat>(x);
  colvec m = as<colvec>(mean);
  int p = xx.n_rows;
  int n = xx.n_cols;
  NumericVector y(n);

  mat Rinv = inv(trimatu(chol(as<arma::mat>(varcovM))));
  //mat R = chol(as<arma::mat>(varcovM));
  double logSqrtDetvarcovM = sum(log(Rinv.diag()));
  double constant = - p*log2pi2;

  for (int i=0; i < n; i++) {
    colvec x_i = xx.col(i) - m;
    rowvec xRinv = trans(x_i)*Rinv;
    //vec xRinv = solve(trimatl(R.t()), x_i);
    double quadform = sum(xRinv%xRinv);
    if (!Log) {
      y(i) = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
    } else{
      y(i) = -0.5*quadform + logSqrtDetvarcovM + constant;
    }
  }

  return y;

}

