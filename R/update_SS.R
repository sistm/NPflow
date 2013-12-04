#Return updated sufficient statistics S of Normal Wishart distribution
#with new data matrix z

update_SS <- function(z, S){
  S_up <- S
  mu0 <- S[["mu"]]
  kappa0 <- S[["kappa"]]
  nu0 <- S[["nu"]]
  lambda0 <- S[["lambda"]]
  
  kappa1 <- kappa0 + 1
  nu1 <- nu0 + 1
  mu1 <- kappa0/(kappa0 + 1)*mu0 + 1/(kappa0 + 1)*z
  lambda1 <- lambda0 + kappa0/(kappa0 + 1)*(z - mu0)%*%t(z - mu0)
  
  S_up[["mu"]] <- mu1
  S_up[["kappa"]] <- kappa1
  S_up[["nu"]] <- nu1
  S_up[["lambda"]] <- lambda1
  
  return(S_up)
}