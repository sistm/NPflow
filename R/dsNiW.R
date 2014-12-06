#'
#' structured Normal inverse Wishart (sNiW) probability density function for multiple inputs
#'@param U_xi a list of length n of observed mean vectors, each of dimension p
#'@param U_psi a list of length n of observed skew vectors of dimension p
#'@param U_Sigma a list of length n of observed covariance matrices, each of dimension p x p
#'@param U_xi0 a list of length K of mean vector parameters for sNiW, each of dimension p
#'@param U_psi0 a list of length K of mean vector parameters for sNiW, each of dimension p
#'@param U_B0 a list of length K of sturcturing matrix parameters for sNiW, each of dimension 2 x 2
#'@param U_Sigma0 a list of length K of covariance matrix parameters for sNiW, each of dimension p x p
#'@param U_df0 a list of length K of degrees of freedom parameters for sNiW, each of dimension p x p
#'@export

mmsNiWlogpdf <- function(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){

    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-(nu0 + d +1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -log(gamma_mv(nu0/2, p=d)) 
                 -1/2*log(det(kronecker(solve(B0), Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(B0, Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }
    
    ml <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
        resml <- mapply(FUN = loglik, 
                        xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                        Lambda0 = U_Sigma0, nu0 = U_df0,
                        MoreArgs=list(xi, psi, Sigma))
#         constnorm <- sum(resml)
#         if(constnorm>0){
#             resml <- resml/constnorm
#         }
        return(resml)
    }
    
    mml <- function(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
        resmml <- mapply(FUN = ml, 
                        xi = U_xi, psi = U_psi, Sigma = U_Sigma,
                        MoreArgs=list(U_xi0, U_psi0, U_B0, U_Sigma0, U_df0))
    }
    
    res <- mml(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0)
    return(res)
    
}

msNiWlogpdf <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
    
    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-(nu0 + d +1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -log(gamma_mv(nu0/2, p=d)) 
                 -1/2*log(det(kronecker(solve(B0), Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(B0, Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }
    

        res <- mapply(FUN = loglik, 
                        xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                        Lambda0 = U_Sigma0, nu0 = U_df0,
                        MoreArgs=list(xi, psi, Sigma))
    
    res <- ml(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0)
    return(res)
    
}

sNiWlogpdf <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
    
    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-(nu0 + d +1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -log(gamma_mv(nu0/2, p=d)) 
                 -1/2*log(det(kronecker(solve(B0), Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(B0, Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }
    
    res <- loglik(xi, psi, Sigma, xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                  Lambda0 = U_Sigma0, nu0 = U_df0
    )
    
    return(res)
    
}

gamma_mv <- function(x,p){
    pi^(p*(p-1)/4)*prod(gamma(x+(1-1:p)/2))
}

digamma_mv <- function(x,p){
    sum(digamma(x+(1-1:p)/2))
}