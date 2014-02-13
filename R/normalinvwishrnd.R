#Sample from a normal inverse Wishart distribution 
#whose parameter are given by the structure hyper
normalinvwishrnd <- function(hyper){
  
  mu0 = hyper[["mu"]]
  kappa0 = hyper[["kappa"]]
  nu0 = hyper[["nu"]]
  lambda0 = hyper[["lambda"]]
  
  # Sample S from an inverse Wishart distribution
  S = invwishrnd(n = nu0, lambda = lambda0)
  
  # Sample mu from a normal distribution
 
  muSupp <- chol(S/kappa0)%*%matrix(rnorm(length(mu0)), 
                                 nrow=length(mu0), ncol=1)
  mu = mu0 + as.vector(muSupp)
  
  
  return(list("S"=S, "mu"=mu))
}
