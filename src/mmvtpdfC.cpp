#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param mean mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param varcovM list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param df vector of length \code{K} of degree of freedom parameters.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension \code{K x n}.
//'
//'@author Boris Hejblum
//'
//'@export
//'@examples
//'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
//'mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000, Log=FALSE)
//'mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10000000, Log=FALSE)
//'
//'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
//'mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000)
//'mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10000000)
//'
//'mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10)
//'mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10)
//'
//'
//'if(require(microbenchmark)){
//' library(microbenchmark)
//' microbenchmark(mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=1, Log=FALSE),
//'                mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)),
//'                         df=c(1), Log=FALSE),
//'                times=10000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
// [[Rcpp::export]]
NumericMatrix mmvtpdfC(const NumericMatrix & x,
                       const NumericMatrix & mean,
                       const List & varcovM,
                       const NumericVector & df,
                       const bool & Log=true){

    mat xx = as<mat>(x);
    mat m = as<mat>(mean);
    int p = xx.n_rows;
    int n = xx.n_cols;
    int K = m.n_cols;
    NumericMatrix y(K,n);

    for(int k=0; k < K; k++){
        mat Rinv = inv(trimatu(chol(as<arma::mat>(varcovM[k]))));
        //mat R = chol(as<arma::mat>(varcovM[k]));
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        colvec mtemp = m.col(k);
        double dftemp = df(k);

        for (int i=0; i < n; i++) {
            colvec x_i = xx.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            //vec xRinv = solve(trimatl(R.t()), x_i);
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2.0) - lgamma(dftemp/2.0) - log(dftemp*M_PI)*p/2.0 ;
            if (!Log) {
                y(k,i) = pow((1.0 + quadform/dftemp),(-(dftemp + p)/2.0))*exp(a+logSqrtDetvarcovM) ;
            } else{
                y(k,i) = (-(dftemp + p)/2.0)*log(1.0 + quadform/dftemp) + a + logSqrtDetvarcovM ;
            }
        }

    }

    return y;

}
