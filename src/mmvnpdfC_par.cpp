#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2;

//' Parallel C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the 
//'data and n the number of data points 
//'@param mean mean vectors matrix of dimension p x K, K being the number of 
//'distributions for which the density probability has to be ealuated
//'@param varcovM list of length K of variance-covariance matrices, 
//'each of dimensions p x p
//'@param logical flag for returning the log of the probability density 
//'function. Defaults is \code{TRUE}
//'@return matrix of densities of dimension K x n
//'@export
//'@examples
//'#NB: change ~/.R/Makevars 'CC' and 'CXX' to set the compiler 
//'#either to 'gcc' and 'g++', or to 'clang' and 'clang++'
//'library(microbenchmark)
//'
//'K=1000
//'d=10
//'z_mat <- NULL
//'m_list <- list()
//'m_mat <- NULL
//'S_list <- list()
//'
//'for(i in 1:K){
//' z_mat <- c(z_mat, c(rep(1.96,d)))
//' m_list <- c(m_list, list(c(rep(-1.5,d))))
//' m_mat <- c(m_mat, c(rep(-1.5,d)))
//' S_list <- c(S_list, list(0.33*diag(d)))
//' }
//' z_mat <- matrix(z_mat, ncol=K, nrow=d)
//' m_mat <- matrix(m_mat, ncol=K, nrow=d)
//'
//'microbenchmark(mvnpdf(x=z_mat, mean=m_list, varcovM=S_list, Log=FALSE),
//'               mmvnpdfC(x=z_mat, mean=m_mat, varcovM=S_list, Log=FALSE),
//'               mmvnpdfC_par(x=z_mat, mean=m_mat, varcovM=S_list, Log=FALSE, 
//'                      ncores = 2),
//'               times=10L)
//'
// [[Rcpp::export]]
NumericMatrix mmvnpdfC_par(NumericMatrix x, 
                       NumericMatrix mean, 
                       List varcovM,
                       bool Log=true,
                       int ncores=1){    
    
    mat xx = as<mat>(x);
    const mat m = as<mat>(mean); 
    const int p = xx.n_rows;
    const int n = xx.n_cols;
    const int K = m.n_cols;
    NumericMatrix y = NumericMatrix(K,n);
    const double constant = - p*log2pi2;

    omp_set_num_threads(ncores);

    for(int k=0; k < K; k++){
        mat S = varcovM[k];
        mat Rinv = inv(trimatu(chol(S)));
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        colvec mtemp = m.col(k);
        

        #pragma omp parallel for shared(mtemp, y, Rinv)
        for (int i=0; i < n; i++) {
            rowvec xRinv = trans(xx.col(i) - mtemp)*Rinv;
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


