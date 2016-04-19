#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2;

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
//'library(microbenchmark)
//'microbenchmark(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'               mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'               mmvnpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), Log=FALSE),
//'               times=1000L)
//'microbenchmark(mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(0, 0),
//'                      varcovM=diag(2), Log=FALSE),
//'               mvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(0, 0),
//'                       varcovM=diag(2), Log=FALSE),
//'               mmvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                        mean=matrix(c(0, 0), nrow=2, ncol=1),
//'                        varcovM=list(diag(2)), Log=FALSE),
//'               times=1000L)
//'microbenchmark(mvnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                      mean=list(c(0,0),c(-1,-1), c(1.5,1.5)),
//'                      varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'               mmvnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                      mean=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
//'                      varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'               times=1000L)
//'
// [[Rcpp::export]]
NumericMatrix mmvnpdfC(NumericMatrix x,
                       NumericMatrix mean,
                       List varcovM,
                       bool Log=true){

    mat xx = as<mat>(x);
    mat m = as<mat>(mean);
    int p = xx.n_rows;
    int n = xx.n_cols;
    int K = m.n_cols;
    NumericMatrix y = NumericMatrix(K,n);
    double constant = - p*log2pi2;

    for(int k=0; k < K; k++){
        mat S = varcovM[k];
        mat Rinv = inv(trimatu(chol(S)));
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        colvec mtemp = m.col(k);

        for (int i=0; i < n; i++) {
            colvec x_i = xx.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
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

