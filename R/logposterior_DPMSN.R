logposterior_DPMSN <- function(z, xi, psi, Sigma, B, hyper, c, m, alpha, n, a, b, diagVar){
    res <- NA
    
    indfull <- which(m!=0)
    mfull <- m[indfull]
    K <- length(indfull)
    if(!is.list(xi)){
        if(is.null(dim(xi))){
            log_vrais <- sum(log(mvsnpdf(x = z, xi = xi, sigma = Sigma, psi = psi)))
            if(!diagVar){
                log_prior_NNiW <-  sum(log(dNNiW(xi, psi, B, Sigma, hyperprior=hyper, log=TRUE)))
            }else{
                log_prior_NNiW <-  0
            }
        } else{
            log_vrais <- sum(log(mvsnpdf(x = z, xi = xi[, c], sigma = Sigma[, , c], psi = psi[, c])))
            log_prior_NNiW <-  sum(dNNiW(xi[,indfull], psi[,indfull], B[,,indfull], Sigma[,,indfull], hyperprior=hyper, log=TRUE))
        }
    }else{
        log_vrais <- sum(log(mvsnpdf(x = z, xi = xi[as.character(c)], 
                                     sigma = Sigma[as.character(c)], psi = psi[as.character(c)])))
        if(!diagVar){
            log_prior_NNiW <-  sum(log(dNNiW(xi[as.character(indfull)], 
                                             psi[as.character(indfull)],
                                             B[as.character(indfull)], 
                                             Sigma[as.character(indfull)], 
                                             hyperprior=hyper, log=TRUE)))
        }else{
            log_prior_NNiW <-  0
        }
    }
    
    log_prior_alpha <- dgamma(alpha, shape=a, scale=1/b, log=TRUE)
    
    log_clustering <- sum(c(lgamma(alpha), K*log(alpha), lgamma(mfull),-lgamma(alpha+n)))
    
    res <- c(log_vrais, log_clustering, log_prior_NNiW, log_prior_alpha)
    
    return(res)
    
}