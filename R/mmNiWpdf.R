
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

  n<-ncol(mu)
  K<-ncol(U_mu0)
  d<-ncol(U_lambda0[[1]])

  res <- matrix(0,nrow=K,ncol=n)

  logpi <- -d/2*log(2*pi)

  for (k in 1:K){

    lambda0_k <- U_lambda0[[k]]
    n0_k <- U_nu0[k]
    kappa0_k <- U_kappa0[[k]]
    mu0_k <- U_mu0[,k]

    # Gamma part
    k_const <- (- n0_k*d/2*log(2)
                + n0_k/2*log(det(lambda0_k))
                - (d*(d-1)/4*log(pi) + sum(lgamma((n0_k+1-1:d)/2)))
                )

    for (i in 1:n){

      S_i <- Sigma[[i]]
      mu_i <- mu[,i]
      x <- mu_i - mu0_k

      logdetS_i <- -(n0_k + d + 1)/2*log(det(S_i))
      logdetSkappa <- -0.5*log(det(S_i/kappa0_k))

      # Exp part
      logexptrace <- -0.5*sum(diag(lambda0_k%*%solve(S_i)))
      quadform <- -0.5*kappa0_k*t(x)%*%solve(S_i)%*%x

      res[k,i] <- logpi + k_const + logdetS_i + logdetSkappa + logexptrace + quadform
    }
  }

  if(!Log){
    res <- exp(res)
  }

  return(res)
}


