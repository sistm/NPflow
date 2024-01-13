#include <RcppArmadillo.h>
#include "lgamma_mvC.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' C++ implementation of multivariate structured Normal inverse Wishart probability density function for multiple inputs
//'
//'@param xi data matrix of dimensions \code{p x n} where columns contain the observed
//'mean vectors.
//'@param psi data matrix of dimensions \code{p x n} where columns contain the observed
//'skew parameter vectors.
//'@param Sigma list of length \code{n} of observed variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param U_xi0 mean vectors matrix of dimension \code{p x K}, \code{K} being the number of
//'distributions for which the density probability has to be evaluated.
//'@param U_psi0 skew parameter vectors matrix of dimension \code{p x K}.
//'@param U_B0 list of length \code{K} of structured scale matrices,
//'each of dimensions \code{p x p}.
//'@param U_Sigma0 list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param U_df0 vector of length \code{K} of degree of freedom parameters.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'@return matrix of densities of dimension \code{K x n}
//'
//'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019) 
//'Sequential Dirichlet Process Mixtures of Multivariate Skew t-distributions for 
//'Model-based Clustering of Flow Cytometry Data. The Annals of Applied Statistics, 
//'13(1): 638-660. <doi: 10.1214/18-AOAS1209>. <arXiv: 1702.04407>. 
//'\url{https://arxiv.org/abs/1702.04407} \doi{10.1214/18-AOAS1209}
//'
//'@export
//'
// [[Rcpp::export]]
NumericMatrix mmsNiWpdfC(const arma::mat & xi,
                         const arma::mat & psi,
                         const List & Sigma,
                         const arma::mat & U_xi0,
                         const arma::mat & U_psi0,
                         const List & U_B0,
                         const List & U_Sigma0,
                         const NumericVector & U_df0,
                         const bool & Log=true){

    int d = xi.n_rows;
    int n = xi.n_cols;
    int K = U_xi0.n_cols;
    NumericMatrix y = NumericMatrix(K,n);

    const double dlog2pi = -d*log(2.0 * M_PI);

    for(int k=0; k < K; k++){
        colvec mu0 = join_vert(U_xi0.col(k), U_psi0.col(k));
        mat lambda0 = U_Sigma0[k];
        mat B0 = U_B0[k];
        mat B0inv = inv_sympd(B0);
        double nu0 = U_df0(k);

        double k_const = (- nu0*d*log(2.0)
                          + nu0*log(det(lambda0))
                          - 2.0*lgamma_mvC(nu0/2.0, d)
        );


        for (int i=0; i < n; i++) {

            colvec mu = join_vert(xi.col(i), psi.col(i)) - mu0;
            mat sigma = Sigma[i];
            mat sigmainv = inv_sympd(sigma);
            mat quadform = mu.t()*kron(B0inv, sigmainv)*mu;

            double logdet_Kron_B0sigma;
            double sign;
            log_det(logdet_Kron_B0sigma, sign, kron(B0, sigma));

            double logdet_sigma;
            double sign2;
            log_det(logdet_sigma, sign2, sigma);

            double res = (-(nu0 + d + 1.0)*logdet_sigma
                          - logdet_Kron_B0sigma
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
