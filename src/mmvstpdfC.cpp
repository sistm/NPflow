#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param xi mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param psi skew parameter vectors matrix of dimension \code{p x K}.
//'@param sigma list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param df vector of length K of degree of freedom parameters.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension \code{K x n}.
//'
//'@author Boris Hejblum
//'
//'@export
//'@examples
//'mmvstpdfC(x = matrix(c(3.399890,-5.936962), ncol=1), xi=matrix(c(0.2528859,-2.4234067)),
//'psi=matrix(c(11.20536,-12.51052), ncol=1),
//'sigma=list(matrix(c(0.2134011, -0.0382573, -0.0382573, 0.2660086), ncol=2)),
//'df=c(7.784106)
//')
//'mvstpdf(x = matrix(c(3.399890,-5.936962), ncol=1), xi=c(0.2528859,-2.4234067),
//'psi=c(11.20536,-12.51052),
//'sigma=matrix(c(0.2134011, -0.0382573, -0.0382573, 0.2660086), ncol=2),
//'df=c(7.784106)
//')
//'
//'#skew-normal limit
//'mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2))
//'          )
//'mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'        xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
//'        df=100000000
//'        )
//'mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)),
//'          df=100000000
//'          )
//'
//'#non-skewed limit
//'mmvtpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'         mean=matrix(c(0, 0)), varcovM=list(diag(2)),
//'         df=10
//'         )
//'mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'          xi=matrix(c(0, 0)), psi=matrix(c(0, 0),ncol=1), sigma=list(diag(2)),
//'          df=10
//'          )
//'
//'if(require(microbenchmark)){
//' library(microbenchmark)
//' microbenchmark(mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                        xi=c(0, 0), psi=c(1, 1),
//'                        sigma=diag(2), df=10),
//'                mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1),
//'                          sigma=list(diag(2)), df=10),
//'                times=1000L)
//'}else{
//' cat("package 'microbenchmark' not available\n")
//'}
// [[Rcpp::export]]
NumericMatrix mmvstpdfC(arma::mat x,
                        arma::mat xi,
                        arma::mat psi,
                        List sigma,
                        NumericVector df,
                        bool Log=true){

    int p = x.n_rows;
    int n = x.n_cols;
    int K = xi.n_cols;
    NumericMatrix y = NumericMatrix(K,n);

    for(int k=0; k < K; k++){
        colvec mtemp = xi.col(k);
        mat psitemp = psi.col(k);
        mat sigmatemp = sigma[k];
        double dftemp = df(k);

        mat omega = sigmatemp + psitemp*trans(psitemp);
        mat omegaInv = inv_sympd(omega);
        mat Rinv=inv(trimatu(chol(omega)));
        mat smallomega = diagmat(sqrt(diagvec(omega)));
        vec alphnum = smallomega*omegaInv*psitemp;
        mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
        vec alphden = rep(alphtemp(0,0), alphnum.size());
        vec alph = alphnum/alphden;
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));

        for (int i=0; i < n; i++) {
            colvec x_i = x.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            mat Qy = x_i.t()*omegaInv*x_i;
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2.0) - lgamma(dftemp/2.0) - log(dftemp*M_PI)*p/2.0;
            double part1 = log(2.0) + (-(dftemp + p)/2.0)*log(1.0 + quadform/dftemp) + a +logSqrtDetvarcovM ;
            //double part1 = 2*pow((1 + quadform/dftemp),(-(dftemp + p)/2))*exp(a+logSqrtDetvarcovM);
            mat quant = trans(alph)*diagmat(1/sqrt(diagvec(omega)))*x_i*sqrt((dftemp + p)/(dftemp+Qy));
            double part2 = ::Rf_pt(quant(0,0), (dftemp + p) , 1, 0);
            if (!Log) {
                y(k,i) = exp(part1)*part2;
            } else{
                y(k,i) = part1 + log(part2);
            }
        }
    }

    return y;

}
