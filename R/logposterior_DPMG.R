#'@keywords internal
#'@importFrom stats dgamma
logposterior_DPMG <- function(z, mu, Sigma, hyper, c, m, alpha, n, a, b){
    res <- NA


    indfull <- which(m!=0)
    mfull <- m[indfull]
    K <- length(indfull)

    if(!is.list(mu)){
        log_vrais <- sum(mvnpdf(x = z, mean = mu[, c], varcovM = Sigma[, , c],Log=TRUE))



        log_prior_NiW <-log(dNiW(mu[,indfull], Sigma[,,indfull], hyperprior=hyper))
    }else{
        log_vrais <- sum(log(mvnpdf(x = z, mean = mu[as.character(c)], varcovM = Sigma[as.character(c)])))
        log_prior_NiW <-  sum(log(dNiW(mu[as.character(indfull)], Sigma[as.character(indfull)], hyperprior=hyper)))
    }

    log_prior_alpha <- stats::dgamma(alpha, shape=a, scale=1/b, log=TRUE)

    log_clustering <- sum(c(lgamma(alpha), K*log(alpha), lgamma(mfull),-lgamma(alpha+n)))

    if(is.list(mu)){

    }else{

    }


    res <- c(log_vrais,log_clustering,log_prior_NiW,log_prior_alpha)

    return(res)

}