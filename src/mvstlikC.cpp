#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' C++ implementation of multivariate skew t likelihood function for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the
//'data and n the number of data points
//'@param c integer vector of cluster allocations with values from 1 to K
//'@param clustval vector of unique values from c in the order corresponding to
//'the storage of cluster parameters in \code{xi}, \code{psi}, and \code{sigma}
//'@param xi mean vectors matrix of dimension p x K, K being the number of
//'clusters
//'@param psi skew parameter vectors matrix of dimension \code{p x K}
//'@param sigma list of length \code{K} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param df vector of length \code{K} of degree of freedom parameters.
//'@param loglik logical flag or returning the log-likelihood intead of the likelihood.
//'Default is \code{TRUE}.
//'@return a list:
//'\itemize{
//'\item{\code{"indiv"}:}{ vector of likelihood of length n;}
//'\item{\code{"clust"}:}{ vector of likelihood of length K;}
//'\item{\code{"total"}:}{ total (log)-likelihood;}
//'}
//'
//'@author Boris Hejblum
//'
// [[Rcpp::export]]
List mvstlikC(arma::mat x,
              arma::vec c,
              arma::vec clustval,
              arma::mat xi,
              arma::mat psi,
              List sigma,
              NumericVector df,
              bool loglik=true){

    int p = x.n_rows;
    int n = x.n_cols;
    int K = xi.n_cols;
    vec yindiv = vec(n);
    vec yclust = vec(K);

    for(int k=0; k < K; k++){
        uvec indk = find(c==clustval(k));
        mat xtemp= x.cols(indk);
        int ntemp = indk.size();

        colvec mtemp = xi.col(k);
        mat psitemp = psi.col(k);
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

        for (int i=0; i < ntemp; i++) {
            colvec x_i = xtemp.col(i) - mtemp;
            rowvec xRinv = trans(x_i)*Rinv;
            mat Qy = x_i.t()*omegaInv*x_i;
            double quadform = sum(xRinv%xRinv);
            double a = lgamma((dftemp + p)/2.0) - lgamma(dftemp/2.0) - log(dftemp*M_PI)*p/2.0;
            double part1 = log(2.0) + (-(dftemp + p)/2)*log(1.0 + quadform/dftemp) + a + logSqrtDetvarcovM;
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
