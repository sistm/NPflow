#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate Normal probability density function for multiple inputs from summarized observations
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param mean mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param varcovM list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param varcov_indiv list of length \code{n} of individual within sum of squares,
//'each of dimensions \code{p x p}.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@param obs_weights a vector for weighting observations in likelihood computation. 
//'@return matrix of densities of dimension \code{K x n}.
//'@export
//'@examples
//'if(require(microbenchmark)){
//'library(microbenchmark)
//' microbenchmark(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
//'                mmvnpdfC(x=matrix(1.96), mean=matrix(0), 
//'                         varcovM=list(diag(1)), Log=FALSE, 
//'                         obs_weights = rep(1, 1)),
//'                times=1000L)
//' microbenchmark(mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
//'                       varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
//'                mvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
//'                        varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
//'                mmvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                         mean=matrix(c(-0.2, 0.3), nrow=2, ncol=1),
//'                         varcovM=list(matrix(c(2, 0.2, 0.2, 2), ncol=2)), 
//'                         Log=FALSE, obs_weights = rep(1, 1)),
//'                times=1000L)
//' microbenchmark(mvnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                       mean=list(c(0,0),c(-1,-1), c(1.5,1.5)),
//'                       varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'                mmvnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                         mean=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
//'                         varcovM=list(diag(2),10*diag(2), 20*diag(2)), 
//'                         Log=FALSE, obs_weights = rep(1, 2)),
//'                times=1000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
// [[Rcpp::export]]
arma::mat mmvnpdfC_summary(arma::mat x,
                           arma::mat mean,
                           List varcov_indiv,
                           List varcovM,
                           arma::rowvec obs_weights,
                           bool Log=true){
    
    int p = x.n_rows;
    int n = x.n_cols;
    int K = mean.n_cols;
    mat y(K,n);
    double constant = - p*log2pi2;
    rowvec wsr = sqrt(obs_weights);
    
    for(int k=0; k < K; k++){
        colvec mtemp = mean.col(k);
        NumericMatrix sigma = varcovM[k];
        arma::mat sigma_arma = arma::mat(sigma.begin(), sigma.nrow(), sigma.ncol(), false);
        
        mat Rinv = inv(trimatu(chol(sigma_arma)));
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        
        for (int i=0; i < n; i++) {
            NumericMatrix vintra = varcov_indiv[i];
            arma::mat vintra_arma = arma::mat(vintra.begin(), vintra.nrow(), vintra.ncol(), false);
            
            colvec x_i = (x.col(i) - mtemp) * wsr(i);
            mat MtM = x_i*trans(x_i);
            mat temp  = (MtM + vintra_arma) % (Rinv*Rinv.t());
            double quadform = accu((MtM + vintra_arma) % (Rinv*Rinv.t())); // sum of all element-wise multiplied elements
            
            y(k,i) = -0.5*quadform + logSqrtDetvarcovM + constant;
        }
    }
    if (!Log) {
        y = exp(y);
    }
    
    return y;
    
}
