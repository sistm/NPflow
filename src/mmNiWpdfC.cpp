#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
double lgamma_mv2C(double x,
                  double p){

    double res = (p*(p-1)/4)*log(M_PI);

    for(int i=0; i < p; i++){
        res = res + lgamma(x - (double)i/2);
    }

    return (res);
}



//' C++ implementation of multivariate Normal inverse Wishart probability density function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the
//'data and n the number of data points
//'@param Mu mean vectors matrix of dimension p x K, K being the number of
//'distributions for which the density probability has to be ealuated
//'@param varcovM list of length K of variance-covariance matrices,
//'each of dimensions p x p
//'@param U_Nu0 vector of length K of degree of freedom parameters
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension K x n
//'@export
//'
//'
// [[Rcpp::export]]
NumericMatrix mmNiWpdfC(NumericMatrix Mu,
                         List Sigma,
                         NumericMatrix U_Mu0,
                         NumericVector U_Kappa0,
                         NumericVector U_Nu0,
                         List U_Sigma0,
                         bool Log=true){

    mat mu = as<mat>(Mu);
    mat mu0 = as<mat>(U_Mu0);

    int d = mu.n_rows;
    int n = mu.n_cols;
    int K = mu0.n_cols;
    NumericMatrix y = NumericMatrix(K,n);

    const double dlog2pi = -d*log(2.0 * M_PI);

    for(int k=0; k < K; k++){

        mat lambda0   = U_Sigma0[k];
        double nu0    = U_Nu0(k);
        double kappa0 = U_Kappa0(k);



        double k_const = (- nu0*d*log(2)
                          + nu0*log(det(lambda0))
                          - 2*lgamma_mv2C(nu0/2, d)
        );


        for (int i=0; i < n; i++) {

            vec mu_ = mu.col(i) - mu0.col(k);
            mat sigma = Sigma[i];
            mat sigmainv = inv_sympd(sigma);
            mat quadform = kappa0*mu_.t()*sigmainv*mu_;

            double logdet_kappa0sigma;
            double sign;
            log_det(logdet_kappa0sigma, sign, sigma/kappa0);

            double logdet_sigma;
            double sign2;
            log_det(logdet_sigma, sign2, sigma);

            double res = (-(nu0 + d +1)*logdet_sigma
                          - logdet_kappa0sigma
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
