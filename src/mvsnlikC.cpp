#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate skew normal likelihood function for multiple inputs
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
//'@param loglik logical flag or returning the log-likelihood instead of the likelihood.
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
List mvsnlikC(const arma::mat & x,
              const arma::vec & c,
              const arma::vec & clustval,
              const arma::mat & xi,
              const arma::mat & psi,
              const List & sigma,
              const bool & loglik=true){

  int p = x.n_rows;
  int n = x.n_cols;
  int K = xi.n_cols;
  vec yindiv = vec(n);
  vec yclust = vec(K);
  double constant = - p*log2pi2;

  for(int k=0; k < K; k++){
    uvec indk = find(c==clustval(k));
    mat xtemp= x.cols(indk);
    int ntemp = indk.size();

    colvec mtemp = xi.col(k);
    mat psitemp = psi.col(k);

    mat omega = as<arma::mat>(sigma[k]) + psitemp*trans(psitemp);
    mat omegaInv =  inv_sympd(omega);
    mat Rinv = inv(trimatu(chol(omega)));
    //mat R = chol(omega);
    mat smallomega = diagmat(sqrt(diagvec(omega)));
    vec alphnum = smallomega*omegaInv*psitemp;
    mat alphtemp =sqrt(1-trans(psitemp)*omegaInv*psitemp);
    vec alphden = rep(alphtemp(0,0), alphnum.size());
    vec alph = alphnum/alphden;
    double logSqrtDetvarcovM = sum(log(Rinv.diag()));

    for (int i=0; i < ntemp; i++) {
      colvec x_i = xtemp.col(i) - mtemp;
      rowvec xRinv = trans(x_i)*Rinv;
      //vec xRinv = solve(trimatl(R.t()), x_i);
      double quadform = sum(xRinv%xRinv);
      double part1 = log(2.0) -0.5*quadform + logSqrtDetvarcovM + constant;
      mat quant = trans(alph)*diagmat(1/sqrt(diagvec(omega)))*x_i;
      double part2 = Rcpp::stats::pnorm_0(quant(0,0), 1, 0);
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
