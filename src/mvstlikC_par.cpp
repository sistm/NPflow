#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
//' Parallel C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the 
//'data and n the number of data points 
//'@param c integer vector of cluster allocations with values from 1 to K
//'@param clustval vector of unique values from c in the order corresponding to 
//'the storage of cluster parameters in \code{xi}, \code{psi}, and \code{varcovM}
//'@param xi mean vectors matrix of dimension p x K, K being the number of 
//'clusters
//'@param psi skew parameter vectors matrix of dimension p x K
//'@param varcovM list of length K of variance-covariance matrices, 
//'each of dimensions p x p
//'@param df vector of length K of degree of freedom parameters
//'@return vector of likelihood of length n
//'@seealso mmvstpdfC, mvstpdf
//'@export
// [[Rcpp::export]]
List mvstlikC_par(NumericMatrix x, 
                  IntegerVector c, 
                  IntegerVector clustval, 
                  NumericMatrix xi, 
                  NumericMatrix psi, 
                  List sigma, 
                  NumericVector df, 
                  bool loglik=true, 
                  int ncores=1){
    
    mat xx = as<mat>(x);
    const mat mxi = as<mat>(xi); 
    const mat mpsi = as<mat>(psi);
    const vec cc = as<vec>(c);
    const vec cclustval = as<vec>(clustval);
    const int p = xx.n_rows;
    const int n = xx.n_cols;
    const int K = mxi.n_cols;
    vec yindiv = vec(n);
    vec yclust = vec(K);
    
    omp_set_num_threads(ncores);
    
    for(int k=0; k < K; k++){
        uvec indk = find(cc==cclustval(k));
        mat xtemp= xx.cols(indk);
        int ntemp = indk.size();
        
        colvec mtemp = mxi.col(k);
        mat psitemp = mpsi.col(k);
        mat sigmatemp = sigma[k];
        double dftemp = df(k);
        
        mat omega = sigmatemp + psitemp*trans(psitemp);
        mat omegaInv =  inv_sympd(omega);
        mat Rinv=inv(trimatu(chol(omega)));
        mat smallomega = diagmat(sqrt(diagvec(omega)));
        vec alphnum = smallomega*omegaInv*psitemp;
        mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
        vec alphden = rep(alphtemp(0,0), alphnum.size());
        vec alph = alphnum/alphden;
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        
        #pragma omp parallel for shared(mtemp, xtemp, dftemp, yindiv, Rinv, omegaInv)
        for (int i=0; i < ntemp; i++) {
            colvec x_i = xtemp.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            mat Qy = x_i.t()*omegaInv*x_i;
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2) - lgamma(dftemp/2) - log(dftemp*M_PI)*p/2;
            double part1 = log(2) + (-(dftemp + p)/2)*log(1 + quadform/dftemp) + a + logSqrtDetvarcovM;
            mat quant = trans(alph)*diagmat(1/sqrt(diagvec(omega)))*x_i*sqrt((dftemp + p)/(dftemp+Qy));
            double part2 = ::Rf_pt(quant(0,0), (dftemp + p) , 1, 0);
            yindiv(indk(i)) = part1 + log(part2);
        }
        yclust(k) = sum(yindiv(indk));
    }
    double ytot = sum(yclust);
    
    if(!loglik){
        yindiv = exp(yindiv);
        yclust = exp(yclust);
        ytot = exp(ytot);
    }
    
    
    return Rcpp::List::create(Rcpp::Named("indiv") = wrap(yindiv),
                              Rcpp::Named("clust") = wrap(yclust),
                              Rcpp::Named("total") = wrap(ytot)
                             );
    
}
