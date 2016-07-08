#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param mean mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param varcovM list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension \code{K x n}.
//'@export
//'@examples
//'if(require(microbenchmark)){
//'library(microbenchmark)
//' microbenchmark(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                mmvnpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), Log=FALSE),
//'                times=1000L)
//' microbenchmark(mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
//'                       varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
//'                mvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
//'                        varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
//'                mmvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                         mean=matrix(c(-0.2, 0.3), nrow=2, ncol=1),
//'                         varcovM=list(matrix(c(2, 0.2, 0.2, 2), ncol=2)), Log=FALSE),
//'                times=1000L)
//' microbenchmark(mvnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                       mean=list(c(0,0),c(-1,-1), c(1.5,1.5)),
//'                       varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'                mmvnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                         mean=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
//'                         varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'                times=1000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
// [[Rcpp::export]]
NumericMatrix mmvnpdfC(arma::mat x,
                       arma::mat mean,
                       List varcovM,
                       bool Log=true){

    int p = x.n_rows;
    int n = x.n_cols;
    int K = mean.n_cols;
    NumericMatrix y(K,n);
    double constant = - p*log2pi2;

    for(int k=0; k < K; k++){
        mat Rinv = inv(trimatu(chol(as<arma::mat>(varcovM[k]))));
        //mat R = chol(as<arma::mat>(varcovM[k]));
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        colvec mtemp = mean.col(k);

        for (int i=0; i < n; i++) {
            colvec x_i = x.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            //vec xRinv = solve(trimatl(R.t()), x_i);
            double quadform = sum(xRinv%xRinv);
            if (!Log) {
                y(k,i) = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
            } else{
                y(k,i) = -0.5*quadform + logSqrtDetvarcovM + constant;
            }
        }
    }

    return y;

}

