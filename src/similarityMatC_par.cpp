#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Parallel C++ implementation
//' 
//'
//'@param c list of MCMC partitions
//'
//'@export
//'
//'@examples
//'c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
//'similarityMatC(sapply(c, "["))
//'
//'c2 <- list()
//'for(i in 1:100){
//'     c2 <- c(c2, list(rmultinom(n=1, size=3000, prob=rexp(n=3000))))
//'}
//'library(microbenchmark)
//'f <- function(){c3 <-sapply(c2, "[")
//'             similarityMatC(c3)}
//'microbenchmark(f(), time=1L)
//'
//'
// [[Rcpp::export]]
List similarityMatC_par(NumericMatrix c,
                        int ncores=1){
    
    const mat cc = as<mat>(c);
    const int N = cc.n_cols;
    const int n = cc.n_rows;
    
    NumericVector cost = NumericVector(N);
    mat similarity = mat(n, n, fill::eye);
    
    omp_set_num_threads(ncores);
    // TODO test several position for #pragma omp parallel for shared(cost)
    
    for(int i=0; i<n-1; i++){ 
        for(int j=i+1; j<n; j++){
            similarity(i,j) = sum(cc.row(i) == cc.row(j))/N;
            similarity(j,i) = similarity(i,j);
            #pragma omp parallel for shared(cost)
            for(int k=0; k<N; k++){
                cost(k) = cost(k) + std::abs(similarity(i,j) - (cc(i,k) == cc(j,k)));
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("similarity") = similarity,
                              Rcpp::Named("cost")=cost*2.0);
}

