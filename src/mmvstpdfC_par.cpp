#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
//' ParallelC++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the 
//'data and n the number of data points 
//'@param xi mean vectors matrix of dimension p x K, K being the number of 
//'distributions for which the density probability has to be ealuated
//'@param psi skew parameter vectors matrix of dimension p x K, K being the number of 
//'distributions for which the density probability has to be ealuated
//'@param varcovM list of length K of variance-covariance matrices, 
//'each of dimensions p x p
//'@param df vector of length K of degree of freedom parameters
//'@param logical flag for returning the log of the probability density 
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension K x n
//'@export
//'@examples
//'
//'library(microbenchmark)
//'
//'K=1000
//'d=10
//'z_mat <- NULL
//'xi_list <- list()
//'xi_mat <- NULL
//'psi_list <- list()
//'psi_mat <- NULL
//'S_list <- list()
//'
//'for(i in 1:K){
//' z_mat <- c(z_mat, c(rep(1.96,d)))
//' xi_list <- c(xi_list, list(c(rep(-1.5,d))))
//' xi_mat <- c(xi_mat, c(rep(-1.5,d)))
//' psi_list <- c(psi_list, list(c(rep(-0.2,d))))
//' psi_mat <- c(psi_mat, c(rep(-0.2,d)))
//' S_list <- c(S_list, list(0.33*diag(d)))
//' }
//' z_mat <- matrix(z_mat, ncol=K, nrow=d)
//' xi_mat <- matrix(xi_mat, ncol=K, nrow=d)
//' psi_mat <- matrix(psi_mat, ncol=K, nrow=d)
//'
//'microbenchmark(mmvstpdfC(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, df=rep(10,K), Log=FALSE),
//'               mmvstpdfC_par(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, df=rep(10,K),  Log=FALSE, 
//'                      ncores = 4),
//'               times=50L)
//'
// [[Rcpp::export]]
NumericMatrix mmvstpdfC_par(NumericMatrix x, 
                        NumericMatrix xi, 
                        NumericMatrix psi, 
                        List sigma, 
                        NumericVector df,
                        bool Log=true, 
                        int ncores=1){
    
    mat xx = as<mat>(x);
    const mat mxi = as<mat>(xi); 
    const mat mpsi = as<mat>(psi); 
    const int p = xx.n_rows;
    const int n = xx.n_cols;
    const int K = mxi.n_cols;
    NumericMatrix y = NumericMatrix(K,n);
    
    omp_set_num_threads(ncores);
    
    for(int k=0; k < K; k++){
        colvec mtemp = mxi.col(k);
        mat psitemp = mpsi.col(k);
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
        
        #pragma omp parallel for shared(mtemp, y, Rinv, omegaInv)
        for (int i=0; i < n; i++) {
            colvec x_i = xx.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            mat Qy = x_i.t()*omegaInv*x_i;
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2) - lgamma(dftemp/2) - log(dftemp*M_PI)*p/2;
            double part1 = log(2) +(-(dftemp + p)/2)*log(1 + quadform/dftemp) + a+logSqrtDetvarcovM ;
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
