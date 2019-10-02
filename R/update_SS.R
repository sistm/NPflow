#'Return updated sufficient statistics S with new data matrix z
#'
#'For internal use only.
#'
#'@param z data matrix
#'
#'@param S previous sufficient statistics
#'
#'@param hyperprior Default is \code{NULL}
#'
#'@keywords internal
#'
#'@importFrom stats weighted.mean
#'
#'@export


update_SS <- function(z, S, hyperprior=NULL, obs_weights = NULL){
  S_up <- S
  mu0 <- S[["mu"]]
  kappa0 <- S[["kappa"]]
  nu0 <- S[["nu"]]
  lambda0 <- S[["lambda"]]
  
  if(length(dim(z))>1 & dim(z)[2]>1){
    
    n <- ncol(z)
    kappa1 <- kappa0 + n
    nu1 <- nu0 + n
    if(is.null(obs_weights)){
      zbar <- apply(X=z, MARGIN=1, FUN=mean)
      unscvarz <- tcrossprod(z-zbar)
    }else{
      zbar <- apply(X = z, MARGIN = 1, FUN = weighted.mean, w = obs_weights)
      unscvarz <- crossprod(t(z-zbar)*sqrt(obs_weights))
    }
    mu1 <- n/(kappa0 + n)*zbar + kappa0/(kappa0 + n)*mu0
    
    # varz <- tcrossprod(z[,1]-zbar)
    # for(j in 2:n){
    #   varz <- varz + tcrossprod(z[,j]-zbar)
    # } 
    # the same computation is achievable though varz <- tcrossprod(z-zbar)
    
    
    if(!is.null(hyperprior)){
      #g0 <- ncol(lambda0) + 5
      g0 <- nu0
      lambda0 <- wishrnd(n = nu0 + g0, 
                         Sigma = solve( solve(lambda0) + solve(hyperprior[["Sigma"]]) )
      )
    }
    lambda1 <- (lambda0 + kappa0*n/(kappa0 + n)*tcrossprod(zbar - mu0) + unscvarz)
    #cat("lambda0 =", lambda0/(nu1-3), "\n")
    #cat("lambda1 =", lambda1/(nu1-3), "\n")
    #cat("varz =", varz/(nu1-3), "\n")
    
  }else{
    kappa1 <- kappa0 + 1
    nu1 <- nu0 + 1
    mu1 <- (kappa0/(kappa0 + 1)*mu0 + 1/(kappa0 + 1)*z)[,1]
    
    if(!is.null(hyperprior)){
      #g0 <- ncol(lambda0) + 5
      g0 <- nu0
      lambda0 <- wishrnd(n = nu0 + g0, 
                         Sigma = solve( solve(lambda0) + solve(hyperprior[["Sigma"]]) )
      )
    }
    lambda1 <- lambda0 + kappa0/(kappa0 + 1)*tcrossprod(z[,1] - mu0)
  }
  
  S_up[["mu"]] <- mu1
  S_up[["kappa"]] <- kappa1
  S_up[["nu"]] <- nu1
  S_up[["lambda"]] <- lambda1
  
  return(S_up)
}