#'@keywords internal
#'@importFrom stats dgamma
logposterior_DPMG <- function(z, mu, Sigma, hyper, c, m, alpha, log_alpha, n, a, b, diagVar){
  
  res <- NA
  
  indfull <- which(m!=0)
  mfull <- m[indfull]
  K <- length(indfull)
  
  
  if(!is.list(mu)){
    U_mu_full <- sapply(indfull, function(j) mu[, j])
    U_Sigma_full <- lapply(indfull, function(j) Sigma[, ,j])
  }else{
    U_mu_full <- sapply(mu, "[")
    U_Sigma_full <- Sigma
  }
  if(nrow(z)==1){
    U_Sigma_full <- lapply(U_Sigma_full, FUN=matrix, nrow=1, ncol=1)
    U_mu_full <- matrix(U_mu_full, nrow=1)
  }
  
  log_lik <- mvnlikC(x=z, c=c, clustval=indfull,
                     mu=U_mu_full, sigma=U_Sigma_full,
                     loglik=TRUE)
  log_vrais <- log_lik$total
  
  
  
  
  if(!diagVar){
    log_prior_NiW <-  sum(dNiW(lapply(indfull, function(j) mu[, j]),
                               U_Sigma_full,
                               hyperprior=hyper, log=TRUE))
  }else{
    betas <- lapply(U_Sigma_full, diag)
    beta0 <- diag(hyper$lambda)
    S <- lapply(betas, function(b){sum(stats::dgamma(x=b,shape=hyper$nu,
                                                     rate=1/beta0, log=TRUE))})
    log_prior_NiW <- sum(unlist(S))
  }
  
  
  lgamma_alpha <- lgamma(alpha)
  
  #cat("truth: ", lgamma_alpha, "\napprox:", lgamma(1 + exp(log_alpha)) - log_alpha, "\n")
  if(is.infinite(lgamma_alpha) | is.nan(lgamma_alpha)){
    lgamma_alpha <- lgamma(1 + exp(log_alpha)) - log_alpha
  }
  log_prior_alpha <- a*log(b) - lgamma(a) + (a-1)*lgamma_alpha - alpha*b
  
  log_clustering <- sum(c(lgamma_alpha, K*log_alpha, lgamma(mfull), -lgamma(alpha+n)))
  
  
  
  res <- c(log_vrais, log_clustering, log_prior_NiW, log_prior_alpha)
  
  return(res)
}


