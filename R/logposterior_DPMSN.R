#'@keywords internal
#'@importFrom stats dgamma
logposterior_DPMSN <- function(z, xi, psi, Sigma, B, hyper, c, m, alpha, n, a, b, diagVar){
    res <- NA

    indfull <- which(m!=0)
    mfull <- m[indfull]
    K <- length(indfull)
    if(!is.list(xi)){
        if(is.null(dim(xi))){
            log_vrais <- sum(mvsnpdf(x = z, xi = xi, sigma = Sigma, psi = psi, Log=TRUE))
            if(!diagVar){
                log_prior_NNiW <-  sum(dNNiW(xi, psi, B, Sigma, hyperprior=hyper, log=TRUE))
            }else{
                log_prior_NNiW <-  0
            }
        }else{
            log_vrais <- sum(mvsnpdf(x = z, xi = xi[, c], sigma = Sigma[, , c], psi = psi[, c], Log=TRUE))

            if(!diagVar){
                log_prior_NNiW <-  sum(dNNiW(xi[,indfull], psi[,indfull], B[,,indfull], Sigma[,,indfull], hyperprior=hyper, log=TRUE))
            }else{
                betas <- apply(X=Sigma[,,indfull], MARGIN=3, diag)
                beta0 <- diag(hyper$lambda)
                S <- apply(betas, MARGIN=2, function(b){sum(stats::dgamma(x=b,shape=hyper$nu,
                                                                   rate=1/beta0, log=TRUE))})
                log_prior_NNiW <- sum(unlist(S))
            }
        }
    }else{
        log_vrais <- sum(log(mvsnpdf(x = z, xi = xi[as.character(c)],
                                     sigma = Sigma[as.character(c)], psi = psi[as.character(c)])))
        if(!diagVar){
            log_prior_NNiW <-  sum(dNNiW(lapply(indfull, function(j) xi[, j]),
                                         lapply(indfull, function(j) psi[, j]),
                                         lapply(indfull, function(j) B[, , j]),
                                         lapply(indfull, function(j) Sigma[, ,j]),
                                         hyperprior=hyper, log=TRUE))
        }else{
            betas <- lapply(Sigma[as.character(indfull)], diag)
            beta0 <- diag(hyper$lambda)
            S <- lapply(betas, function(b){sum(stats::dgamma(x=b,shape=hyper$nu,
                                                      rate=1/beta0, log=TRUE))})
            log_prior_NNiW <- sum(unlist(S))
        }
    }

    log_prior_alpha <- stats::dgamma(alpha, shape=a, scale=1/b, log=TRUE)

    log_clustering <- sum(c(lgamma(alpha), K*log(alpha), lgamma(mfull),-lgamma(alpha+n)))

    res <- c(log_vrais, log_clustering, log_prior_NNiW, log_prior_alpha)

    return(res)

}