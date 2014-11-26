#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2;

//' Parallel C++ implementation of multivariate skew Normal probability density function for multiple inputs
//' 
//'
//'@param x data matrix of dimension p x n, p being the dimension of the 
//'data and n the number of data points 
//'@param xi mean vectors matrix of dimension p x K, K being the number of 
//'distributions for which the density probability has to be ealuated
//'@param psi skew parameter vectors matrix of dimension p x K, K being the number of 
//'distributions for which the density probability has to be ealuated
//'@param varcovM list of length K of variance-covariance matrices, 
//'each of dimensions p x p
//'@param logical flag for returning the log of the probability density 
//'function. Defaults is \code{TRUE}.
//'@param ncores integer number of cores to be used in parallel
//'@return matrix of densities of dimension K x n
//'@export 
//'@examples
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
//'microbenchmark(mvsnpdf(x=z_mat, xi=xi_list, psi=psi_list, sigma=S_list, Log=FALSE),
//'               mmvsnpdfC(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, Log=FALSE),
//'               mmvsnpdfC_par(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, Log=FALSE, 
//'                      ncores = 2),
//'               times=10L)
//'microbenchmark(mmvsnpdfC(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, Log=FALSE),
//'               mmvsnpdfC_par(x=z_mat, xi=xi_mat, psi=psi_mat, sigma=S_list, Log=FALSE, 
//'                      ncores = 4),
//'               times=100L)
//'               
//'              
// [[Rcpp::export]]
NumericMatrix mmvsnpdfC_par(NumericMatrix x, 
                        NumericMatrix xi, 
                        NumericMatrix psi, 
                        List sigma,
                        bool Log=true,
                        int ncores=1){
    
    mat xx = as<mat>(x);
    const mat mxi = as<mat>(xi); 
    const mat mpsi = as<mat>(psi); 
    const int p = xx.n_rows;
    const int n = xx.n_cols;
    const int K = mxi.n_cols;
    NumericMatrix y = NumericMatrix(K,n);
    const double constant = - p*log2pi2;
    
    omp_set_num_threads(ncores);
    
    for(int k=0; k < K; k++){
        colvec mtemp = mxi.col(k);
        mat psitemp = mpsi.col(k);
        mat sigmatemp = sigma[k];
        
        mat omega = sigmatemp + psitemp*trans(psitemp);
        mat omegaInv = inv(omega);
        mat Rinv=inv(trimatu(chol(omega)));
        mat smallomega = diagmat(sqrt(diagvec(omega)));
        vec alphnum = smallomega*omegaInv*psitemp;
        mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
        vec alphden = rep(alphtemp(0,0), alphnum.size());
        vec alph = alphnum/alphden;
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        
        #pragma omp parallel for shared(mtemp, y, Rinv)
        for (int i=0; i < n; i++) {
            colvec x_i = xx.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            double quadform = sum(xRinv%xRinv);
            double part1 = log(2) -0.5*quadform + logSqrtDetvarcovM + constant;
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

