#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of residual trace computation step used when sampling the scale
//'
//'@param eps a numeric matrix where each column contains the centered and unskewed observations
//'@param sigma a numeric covariance matrix
//'
//'@return the computed trace
//'
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
NumericVector traceEpsC(arma::mat eps,
                        arma::mat sigma){

    int n =eps.n_cols;
    vec tra(n);

    for(int i=0; i<n; i++){
        colvec eps_i = eps.col(i);
        tra(i) = trace(eps_i*eps_i.t()*inv_sympd(sigma));
    }

    return(wrap(tra));
}
