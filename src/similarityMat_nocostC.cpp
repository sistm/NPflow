#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation
//'
//'
//'@param cc a matrix whose columns each represents a ()MCMC) partition
//'
//'@export
//'
//'@examples
//'c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
//'similarityMat_nocostC(sapply(c, "["))
//'
//'c2 <- list()
//'for(i in 1:10){
//'     c2 <- c(c2, list(rmultinom(n=1, size=1000, prob=rexp(n=1000))))
//'}
//'
//'c3 <- sapply(c2, "[")
//'library(microbenchmark)
//'microbenchmark(similarityMat(c3), similarityMat_nocostC(c3), times=2L)
//'
// [[Rcpp::export]]
List similarityMat_nocostC(arma::mat cc){

    int N = cc.n_cols;
    int n = cc.n_rows;

    NumericVector cost = NumericVector(N);
    mat similarity = mat(n, n, fill::eye);

    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){
            similarity(i,j) = sum(cc.row(i) == cc.row(j))/N;
            similarity(j,i) = similarity(i,j);
        }
    }
    return Rcpp::List::create(Rcpp::Named("similarity") = similarity);
}
