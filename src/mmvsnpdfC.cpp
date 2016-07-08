#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate skew Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param xi mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param psi skew parameter vectors matrix of dimension \code{p x K}.
//'@param sigma list of length K of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension \code{K x n}.
//'
//'@author Boris Hejblum
//'
//'@export
//'@examples
//'mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)), Log=FALSE
//'          )
//'mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2))
//'          )
//'
//'if(require(microbenchmark)){
//' library(microbenchmark)
//' microbenchmark(mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), xi=c(0, 0), psi=c(1, 1),
//'                        sigma=diag(2), Log=FALSE),
//'                mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), xi=matrix(c(0, 0)),
//'                          psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)), Log=FALSE),
//'                times=1000L
//'              )
//' microbenchmark(mvsnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                       xi=list(c(0,0),c(-1,-1), c(1.5,1.5)),
//'                       psi=list(c(0.1,0.1),c(-0.1,-1), c(0.5,-1.5)),
//'                       sigma=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'                mmvsnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
//'                          xi=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
//'                          psi=matrix(c(0.1,0.1,-0.1,-1, 0.5,-1.5), nrow=2, ncol=3),
//'                          sigma=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
//'               times=1000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
// [[Rcpp::export]]
NumericMatrix mmvsnpdfC(arma::mat x,
                        arma::mat xi,
                        arma::mat psi,
                        List sigma,
                        bool Log=true){


    int p = x.n_rows;
    int n = x.n_cols;
    int K = xi.n_cols;
    NumericMatrix y(K,n);
    double constant = - p*log2pi2;

    for(int k=0; k < K; k++){
        colvec mtemp = xi.col(k);
        mat psitemp = psi.col(k);

        mat omega = as<arma::mat>(sigma[k]) + psitemp*trans(psitemp);
        mat omegaInv = inv(omega);
        mat Rinv=inv(trimatu(chol(omega)));
        //mat R = chol(omega);
        mat smallomega = diagmat(sqrt(diagvec(omega)));
        vec alphnum = smallomega*omegaInv*psitemp;
        mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
        vec alphden = rep(alphtemp(0,0), alphnum.size());
        vec alph = alphnum/alphden;
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));

        for (int i=0; i < n; i++) {
            colvec x_i = x.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            //vec xRinv = solve(trimatl(R.t()), x_i);
            double quadform = sum(xRinv%xRinv);
            double part1 = log(2.0) -0.5*quadform + logSqrtDetvarcovM + constant;
            mat quant = trans(alph)*diagmat(1/sqrt(diagvec(omega)))*x_i;
            double part2 = Rcpp::stats::pnorm_0(quant(0,0), 1, 0);
            if (!Log) {
                y(k,i) = exp(part1)*part2;
            } else{
                y(k,i) = part1 + log(part2);
            }
        }
    }

    return y;

}

