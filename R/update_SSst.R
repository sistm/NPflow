#' Return updated sufficient statistics S of Normal Wishart distribution
#' with new data matrix z
#' 
#' For internal use only.
#' 
#'@keywords internal
#'
#'@export

# w 

update_SSst <- function(z, S, ltn, scale, df){
    
    b0_xi <- S[["b_xi"]]
    b0_psi <- S[["b_psi"]]
    D0_xi <- S[["D_xi"]]
    D0_psi <- S[["D_psi"]]
    nu0 <- S[["nu"]]
    lambda0 <- S[["lambda"]]
    
    sc_sr <- sqrt(scale)
    
    if(is.null(dim(z))){
        z <- matrix(z, ncol=1)
    }
    n <- ncol(z)

    X <- matrix(c(sc_sr, sc_sr*ltn), ncol=2, byrow=FALSE)
    B <- solve(crossprod(X)+diag(c(1/D0_xi, 1/D0_psi)))
    temp <- apply(X=z, MARGIN=1, FUN=function(x){x*sc_sr})
    if(n<2){
        temp <- matrix(temp, ncol=nrow(z))
    }
    b <- (crossprod(temp,X) + cbind(b0_xi/D0_xi, b0_psi/D0_psi))%*%B

    
    b_xi <- b[,1]
    b_psi <- b[,2]
    
    nu1 <- nu0 + n #c
    
    eps2 <- tcrossprod(z[,1] - b_xi - ltn[1]*b_psi)*scale[1]
    
    if(n>1){
        for (i in 2:n){
            eps2 <- eps2 + tcrossprod(z[,i] - b_xi - ltn[i]*b_psi)*scale[i]
        }
    }
    
    #conjugate hyperprior on lambda: whishart distribution
    g0 <- ncol(lambda0)
    lambda0 <- wishrnd(n=g0, Sigma=solve(lambda0))
    
    lambda1 <- lambda0 + (eps2 + tcrossprod(b_xi-b0_xi)/D0_xi + tcrossprod(b_psi-b0_psi)/D0_psi)
    
    
    S_up <- list()
    S_up[["b_xi"]] <- b_xi
    S_up[["b_psi"]] <- b_psi
    S_up[["B"]] <- B
    S_up[["nu"]] <- nu1
    S_up[["lambda"]] <- lambda1
    S_up[["D_xi"]] <- D0_xi
    S_up[["D_psi"]] <- D0_psi
    S_up[["df"]] <- df
    
    
    return(S_up)
}