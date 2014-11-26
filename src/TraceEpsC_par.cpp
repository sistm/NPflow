#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Parallel C++ implementation of residual trace computation step used when sampling the scale
//' 
//'@param eps
//'@keywords internal
//'
//'@export
//'
// [[Rcpp::export]]
NumericVector traceEpsC_par(NumericMatrix eps, 
                        NumericMatrix sigma,
                        int ncores=1){
    
    const mat epsm =as<mat>(eps);
    const mat S = as<mat>(sigma);
    const int n =epsm.n_cols;
    vec tra = vec(n);
    
    omp_set_num_threads(ncores);
    
    #pragma omp parallel for shared(tra)
    for(int i=0; i<n; i++){
        colvec eps_i = epsm.col(i);
        tra(i) = trace(eps_i*eps_i.t()*inv_sympd(S));
    }
    
    return(wrap(tra));
}
