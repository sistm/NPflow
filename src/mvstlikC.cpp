#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the 
//'data and n the number of data points 
//'@param c integer vector of cluster allocations with values from 1 to K
//'@param clustval vector unique values of c in the order corresponding to 
//'the storage of cluster parameters in \code{xi}, \code{psi}, and \code{varcovM}
//'@param xi mean vectors matrix of dimension p x K, K being the number of 
//'clusters
//'@param psi skew parameter vectors matrix of dimension p x K
//'@param varcovM list of length K of variance-covariance matrices, 
//'each of dimensions p x p
//'@param df vector of length K of degree of freedom parameters
//'@return vector of likelihood of length n
//'@export
//'@examples
//'
//'#skew-normal limit
//'mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), 
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2))
//'          )
//'mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'        xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
//'        df=100000000
//'        )
//'mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), 
//'          xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)),
//'          df=100000000
//'          )
//'          
//'#non-skewed limit         
//'mmvtpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'         mean=matrix(c(0, 0)), varcovM=list(diag(2)),
//'         df=10
//'         )
//'mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), 
//'          xi=matrix(c(0, 0)), psi=matrix(c(0, 0),ncol=1), sigma=list(diag(2)),
//'          df=10
//'          )
//'          
//'microbenchmark(mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), 
//'                       xi=c(0, 0), psi=c(1, 1), 
//'                       sigma=diag(2), df=10),
//'               mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
//'                         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), 
//'                         sigma=list(diag(2)), df=10),
//'               times=10000L)
//'
// [[Rcpp::export]]
List mvstlikC(NumericMatrix x, IntegerVector c, IntegerVector clustval, NumericMatrix xi, 
NumericMatrix psi, List sigma, NumericVector df, bool loglik){
    
    mat xx = as<mat>(x);
    mat mxi = as<mat>(xi); 
    mat mpsi = as<mat>(psi);
    vec cc = as<vec>(c);
    vec cclustval = as<vec>(clustval);
    int p = xx.n_rows;
    int n = xx.n_cols;
    int K = mxi.n_cols;
    vec yindiv = vec(n);
    vec yclust = vec(K);
    
    for(int k=0; k < K; k++){
        uvec indk = find(cc==cclustval(k));
        mat xtemp= xx.cols(indk);
        int ntemp = indk.size();
        
        colvec mtemp = mxi.col(k);
        mat psitemp = mpsi.col(k);
        mat sigmatemp = sigma[k];
        double dftemp = df(k);
        
        mat omega = sigmatemp + psitemp*trans(psitemp);
        mat omegaInv = inv(omega);
        mat Rinv=inv(trimatu(chol(omega)));
        mat smallomega = diagmat(sqrt(diagvec(omega)));
        vec alphnum = smallomega*omegaInv*psitemp;
        mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
        vec alphden = rep(alphtemp(0,0), alphnum.size());
        vec alph = alphnum/alphden;
        double logSqrtDetvarcovM = sum(log(Rinv.diag()));
        
        for (int i=0; i < ntemp; i++) {
            colvec x_i = xtemp.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            mat Qy = x_i.t()*omegaInv*x_i;
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2) - lgamma(dftemp/2) - log(dftemp*M_PI)*p/2;
            double part1 = 2*pow((1 + quadform/dftemp),(-(dftemp + p)/2))*exp(a+logSqrtDetvarcovM);
            mat quant = trans(alph)*diagmat(1/sqrt(diagvec(omega)))*x_i*sqrt((dftemp + p)/(dftemp+Qy));
            double part2 = ::Rf_pt(quant(0,0), (dftemp + p) , 1, 0);
            yindiv(indk(i)) = log(part1*part2);
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
