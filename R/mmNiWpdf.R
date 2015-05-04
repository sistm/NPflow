
#'multivariate Normal inverse Wishart probability density function for multiple inputs

#'@param x data matrix of dimension p x n, p being the dimension of the 
#'data and n the number of data points 
#'@param Mu mean vectors matrix of dimension p x K, K being the number of 
#'distributions for which the density probability has to be ealuated
#'@param varcovM list of length K of variance-covariance matrices, 
#'each of dimensions p x p
#'@param U_Nu0 vector of length K of degree of freedom parameters
#'@param logical flag for returning the log of the probability density 
#'function. Defaults is \code{TRUE}.
#'@return matrix of densities of dimension K x n
#'@export


mmNiWpdf<-function(mu,Sigma,U_mu0,U_kappa0,U_nu0,U_lambda0,Log){
  K<-ncol(mu)
  J<-ncol(U_mu0)
  res<-matrix(0,nrow=J,ncol=K)
  d<-ncol(U_lambda0[[1]])
  for (k in 1:K){
    S_k<-Sigma[[k]]
    mu_k<-mu[,k]
    for (j in 1:J){
      n0_j<-U_nu0[j]
      lambda0_j<-U_lambda0[[j]]
      kappa0_j<-U_kappa0[[j]]
      mu0_j<-U_mu0[,j]
      
      logpi<- -d/2*log(2*pi)
      logdetS<- -(n0_j+d+1)/2*log(det(S_k))
      
      # Gamma part
      plst<-1:d
      gamln<-lgamma((n0_j+1-plst)/2)
      sumgamln<-sum(gamln)
      logtwo<- -n0_j*d/2*log(2)
      logdetlambda<-n0_j/2*log(det(lambda0_j))
      gampart<- logtwo+logdetlambda-sumgamln
      
      logdetSkappa<- -.5*log(det(S_k/kappa0_j))
      
      # Exp part
      logexptrace<--.5*sum(diag(lambda0_j%*%solve(S_k)))
      x<-mu_k-mu0_j
      quadform<- -kappa0_j/2*t(x)%*%solve(S_k)%*%x
      
      res[j,k]<-logpi+logdetS+gampart+logdetSkappa+logexptrace+quadform
    }
  }
  if(Log){
    return(res)
  }
  else{
    return(exp(res))
  }
}


