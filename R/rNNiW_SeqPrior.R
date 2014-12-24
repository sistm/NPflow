#' Sample from a normal inverse Wishart distribution 
#' whose parameter are given by the structure SufStat
#'
#' 
#' For internal use only.
#' 
#'@keywords internal
#'
#'@export
#'
rNNiW_SeqPrior <- function(SufStat, diagVar){
    b0_xi = SufStat[["b_xi"]]
    b0_psi = SufStat[["b_psi"]]
    B = SufStat[["B"]] # B
    nu0 = SufStat[["nu"]] #c
    lambda0 = SufStat[["lambda"]] #C
    
    if(is.null(B)){
        B <- diag(c(1/SufStat[["D_xi"]], 1/SufStat[["D_psi"]]))
    }else{
        B <- solve(B)
    }
    
    # Sample S from an inverse Wishart distribution
    
    if(diagVar){
        betas <- diag(lambda0)
        S <- diag(1/rgamma(n=length(b0_xi), shape=nu0, 
                           rate=betas))
    }else{
        S = invwishrnd(n = nu0, lambda = lambda0)
    }
    
    # Sample mu from a normal distribution
    d <- length(b0_xi)
    muSupp <- (c(b0_xi, b0_psi) 
               + matrix(rnorm(2*d), nrow=1, ncol=2*d)%*%chol(B%x%S))
    xi <- muSupp[1:d]
    psi <- muSupp[(d+1):(2*d)]
    
    
    return(list("S"=S, "xi"=xi, "psi"=psi))
}

