#Return updated sufficient statistics S of Normal Wishart distribution
#with new data matrix z

update_SS <- function(z, S){
  S_up <- S
  mu0 <- S[["mu"]]
  kappa0 <- S[["kappa"]]
  nu0 <- S[["nu"]]
  lambda0 <- S[["lambda"]]
    
  if(!is.null(dim(z))){
      n <- ncol(z)
      zbar <- apply(X=z, MARGIN=1, FUN=mean)
      
      kappa1 <- kappa0 + n
      nu1 <- nu0 + n
      mu1 <- n/(kappa0 + n)*zbar + kappa0/(kappa0 + n)*mu0
      varz <- (z[,1]-zbar)%*%t(z[,1]-zbar)
      for(j in 2:n){
          varz <- varz + (z[,j]-zbar)%*%t(z[,j]-zbar)
      }
      lambda1 <- (lambda0 + kappa0*n/(kappa0 + n)*(zbar - mu0)%*%t(zbar - mu0)
                  + varz)
  } else{
      kappa1 <- kappa0 + 1
      nu1 <- nu0 + 1
      mu1 <- kappa0/(kappa0 + 1)*mu0 + 1/(kappa0 + 1)*z
      lambda1 <- lambda0 + kappa0/(kappa0 + 1)*(z - mu0)%*%t(z - mu0)   
  }
  
  S_up[["mu"]] <- mu1
  S_up[["kappa"]] <- kappa1
  S_up[["nu"]] <- nu1
  S_up[["lambda"]] <- lambda1
  
  return(S_up)
}