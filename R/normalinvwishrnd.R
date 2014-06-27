#' Sample from a normal inverse-Wishart distribution 
#' whose parameter are given by the structure hyper
#'
#' 
#' For internal use only.
#' 
#'@keywords internal
#'
#'@export normalinvwishrnd
#'

normalinvwishrnd <- function(hyper){
  
  mu0 = hyper[["mu"]]
  kappa0 = hyper[["kappa"]]
  nu0 = hyper[["nu"]]
  lambda0 = hyper[["lambda"]]
  
  # Sample S from an inverse Wishart distribution
  S = invwishrnd(n = nu0, lambda = lambda0)
  
  # Sample mu from a normal distribution
 
  muSupp <- matrix(rnorm(length(mu0)), nrow=1, 
                   ncol=length(mu0))%*%chol(S/kappa0)
  mu = mu0 + as.vector(muSupp)
  
  
  return(list("S"=S, "mu"=mu))
}
