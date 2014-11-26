#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation
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
// [[Rcpp::export]]
List similarityMatC(NumericMatrix c){
    
    mat cc = as<mat>(c);
    int N = cc.n_cols;
    int n = cc.n_rows;
    
    NumericVector cost = NumericVector(N);
    mat similarity = mat(n, n, fill::eye);
    
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            similarity(i,j) = sum(cc.row(i) == cc.row(j))/N;
            similarity(j,i) = similarity(i,j);
            for(int k=0; k<N; k++){
                cost(k) = cost(k) + std::abs(similarity(i,j) - (cc(i,k) == cc(j,k)));
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("similarity") = similarity,
                              Rcpp::Named("cost")=cost*2.0);
}
