#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate log gamma function
//' 
//'@param x strictly positive real number
//'@param p integer
//'
//'@export
//'
// [[Rcpp::export]]
double lgamma_mvC(double x, 
                  double p){
                         
    double res = (p*(p-1)/4)*log(M_PI);
    
    for(int i=0; i < p; i++){
        res = res + lgamma(x - (double)i/2);
    }
    
    return (res);
}



//' C++ implementation of multivariate structured Normal inverse Wishart probability density function for multiple inputs
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
//'
//'
// [[Rcpp::export]]
NumericMatrix mmsNiWpdfC(NumericMatrix xi, 
                         NumericMatrix psi, 
                         List Sigma,
                         NumericMatrix U_xi0, 
                         NumericMatrix U_psi0, 
                         List U_B0, 
                         List U_Sigma0,
                         NumericVector U_df0,
                         bool Log=true){
    
    mat mxi = as<mat>(xi); 
    mat mpsi = as<mat>(psi);
    mat mxi0 = as<mat>(U_xi0); 
    mat mpsi0 = as<mat>(U_psi0); 
    int d = mxi.n_rows;
    int n = mxi.n_cols;
    int K = mxi0.n_cols;
    NumericMatrix y = NumericMatrix(K,n);
    
    const double dlog2pi = -d*log(2.0 * M_PI);
    
    for(int k=0; k < K; k++){
        colvec mu0 = join_vert(mxi0.col(k), mpsi0.col(k));
        mat lambda0 = U_Sigma0[k];
        mat B0 = U_B0[k];
        mat B0inv = inv_sympd(B0);
        double nu0 = U_df0(k);
        
        double k_const = (- nu0*d*log(2) 
                          + nu0*log(det(lambda0)) 
                          - 2*lgamma_mvC(nu0/2, d)
        );
        
        
        for (int i=0; i < n; i++) {
            
            colvec mu = join_vert(mxi.col(i), mpsi.col(i)) - mu0;
            mat sigma = Sigma[i];
            mat sigmainv = inv_sympd(sigma);
            mat quadform = mu.t()*kron(B0inv, sigmainv)*mu;
            
            double logdet_Kron_B0sigma;
            double sign;
            log_det(logdet_Kron_B0sigma, sign, kron(B0, sigma));
            
            double logdet_sigma;
            double sign2;
            log_det(logdet_sigma, sign2, sigma);
            
            double res = (-(nu0 + d +1)*logdet_sigma
                          - logdet_Kron_B0sigma
                          - trace(lambda0*sigmainv)
                          - quadform(0,0)
            );
            
            y(k,i) = (dlog2pi + k_const + res)/2;
            
            if (!Log) {
                y(k,i) = exp(y(k,i));
            }
        }
    }
    
    return y;
    
}
