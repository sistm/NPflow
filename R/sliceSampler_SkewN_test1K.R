sliceSampler_SkewN_test1K <- function(c, m, alpha, z, hyperG0, U_xi, U_psi, U_Sigma){
    
    # latent truncated normal variables
    siginv <- solve(U_Sigma)
    psi <- U_psi
    A_k <-  1/(1 + (crossprod(psi, siginv)%*%psi))
    a_ik <- (tcrossprod(A_k, psi)%*%siginv%*%(z-U_xi))
    ltn <- rtruncnorm(ncol(z), a=0, b=Inf, mean = a_ik, sd = sqrt(A_k))
    
    
    return(list("latentTrunc"=ltn))
}