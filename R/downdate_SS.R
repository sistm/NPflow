#Return updated sufficient statistics S of Normal Wishart distribution 
#by removing data point z

downdate_SS <- function(z, S){
  S_up <- S
  mu1 <- S[["mu"]]
  kappa1 <- S[["kappa"]]
  nu1 <- S[["nu"]]
  lambda1 <- S[["lambda"]]
  
  kappa0 <- kappa1 - 1
  nu0 <- nu1 - 1
  mu0 <- (kappa0 + 1)/kappa0*mu1 - 1/kappa0*z
  lambda0 <- lambda1 - kappa0/(kappa0 + 1)*(z - mu0)%*%t(z - mu0)
  
  S_up[["mu"]] <- mu0
  S_up[["kappa"]] <- kappa0
  S_up[["nu"]] <- nu0
  S_up[["lambda"]] <- lambda0
  
  return(S_up)
}