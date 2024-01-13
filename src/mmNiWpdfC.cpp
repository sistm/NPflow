#include <RcppArmadillo.h>
#include "lgamma_mvC.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate Normal inverse Wishart probability density function for multiple inputs
//'
//'@param Mu data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points, where each column is an observed mean vector.
//'@param Sigma list of length \code{n} of observed variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param U_Mu0 mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated
//'@param U_Sigma0 list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param U_Kappa0 vector of length \code{K} of scale parameters.
//'@param U_Nu0 vector of length \code{K} of degree of freedom parameters.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension K x n
//'
//'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019) 
//'Sequential Dirichlet Process Mixtures of Multivariate Skew t-distributions for 
//'Model-based Clustering of Flow Cytometry Data. The Annals of Applied Statistics, 
//'13(1): 638-660. <doi: 10.1214/18-AOAS1209>. <arXiv: 1702.04407>. 
//'\url{https://arxiv.org/abs/1702.04407} \doi{10.1214/18-AOAS1209}
//'
//'@export
//'
//'
// [[Rcpp::export]]
NumericMatrix mmNiWpdfC(const arma::mat & Mu,
                        const List & Sigma,
                        const arma::mat & U_Mu0,
                        const NumericVector & U_Kappa0,
                        const NumericVector & U_Nu0,
                        const List & U_Sigma0,
                        const bool & Log=true){

  int d = Mu.n_rows;
  int n = Mu.n_cols;
  int K = U_Mu0.n_cols;
  NumericMatrix y(K,n);

  const double dlog2pi = -d*log(2.0 * M_PI);

  for(int k=0; k < K; k++){

    mat lambda0   = U_Sigma0[k];
    double nu0    = U_Nu0(k);
    double kappa0 = U_Kappa0(k);



    double k_const = (- nu0*d*log(2.0)
                        + nu0*log(det(lambda0))
                        - 2.0*lgamma_mvC(nu0/2.0, d)
    );


    for (int i=0; i < n; i++) {

      vec mu_ = Mu.col(i) - U_Mu0.col(k);
      mat sigma = Sigma[i];
      mat sigmainv = inv_sympd(sigma);
      mat quadform = kappa0*mu_.t()*sigmainv*mu_;

      double logdet_kappa0sigma;
      double sign;
      log_det(logdet_kappa0sigma, sign, sigma/kappa0);

      double logdet_sigma;
      double sign2;
      log_det(logdet_sigma, sign2, sigma);

      double res = (-(nu0 + d + 1.0)*logdet_sigma
                      - logdet_kappa0sigma
                      - trace(lambda0*sigmainv)
                      - quadform(0,0)
      );

      y(k,i) = (dlog2pi + k_const + res)/2.0;

      if (!Log) {
        y(k,i) = exp(y(k,i));
      }
    }
  }

  return y;
}
