#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of residual trace computation step used when sampling the scale
//' 
//'@param eps
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
NumericVector traceEpsC(NumericMatrix eps, 
                        NumericMatrix sigma){
    
    mat epsm =as<mat>(eps);
    mat S = as<mat>(sigma);
    int n =epsm.n_cols;
    vec tra = vec(n);
    
    for(int i=0; i<n; i++){
        colvec eps_i = epsm.col(i);
        tra(i) = trace(eps_i*eps_i.t()*inv_sympd(S));
    }
    
    return(wrap(tra));
}
