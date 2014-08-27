sample_scale_test1K <- function(c, m, z, U_xi, U_psi, 
                         U_Sigma, U_df, ltn, weights, scale){
    
    n <- length(c)
    
    
                prior_df <- function(x, log=FALSE){
                    if(log){
                        y <- log(x/100)-x/10
                    }else{
                        y <- x/100*exp(-x/10)
                    }
                    return(y)
                }
    
    #     #Hierarchical prior
    #     prior_df <- function(x, d, log=FALSE){
    #         if(log){
    #             y <- log(2)+log(d)+log(x)-3*log(x+d)
    #         }else{
    #             y <- 2*d*x/(x+d)^3
    #         }
    #         return(y)
    #     }
    
    #Jeffrey's prior
#     prior_df <- function(x, log=FALSE){
#         smalltemp <- trigamma(x/2)-trigamma((x+1)/2)-2*(x+3)/(x*(x+1)^2)
#         if(smalltemp <= 0){
#             y=ifelse(log,-Inf,0)
#         }else{
#             if(log){
#                 y <- 1/2*(log(x)-log(x+3) +log(smalltemp))
#             }else{
#                 y <- sqrt(x/(x+3)*smalltemp)
#             }
#         }
#         return(y)
#     }
    
    
    c_df <- 2
    acc_rate <- 0
    
    df <- U_df
    df_new <- 1+exp(runif(1, min=log(df-1)-c_df, max = log(df-1)+c_df))
    
    loglikold <- mvstlikC(x=z,c=c, clustval=0, xi=matrix(ncol=1,U_xi), psi=matrix(ncol=1,U_psi), sigma =list(U_Sigma), df=df, loglik=TRUE)$total
    logliknew <- mvstlikC(x=z,c=c, clustval=0, xi=matrix(ncol=1,U_xi), psi=matrix(ncol=1,U_psi), sigma =list(U_Sigma), df=df_new, loglik=TRUE)$total
    u <- runif(1)
    
    prob_new <- exp(logliknew + prior_df(df_new, log=TRUE) + log(df_new-1)
                    -(loglikold + prior_df(df, log=TRUE) + log(df-1)))
    if(prob_new==-Inf){
        prob_new <- 0
    }else if(prob_new==Inf){
        prob_new <- 1
    }else{ 
        prob_new <- min(1, prob_new)
    }
    if (u<prob_new){
        acc_rate <- acc_rate + 1
        df <- df_new
    }
    
    
    eps <- z - U_xi - sapply(X=ltn, FUN=function(x){x*U_psi})
    tra <- apply(X=eps, MARGIN=2, FUN=function(v){sum(diag(tcrossprod(v)%*%solve(U_Sigma)))})
    
    scale <- rgamma(n, shape=(df + nrow(z) + 1)/2, 
                    rate=(df + ltn^2 + tra)/2)
    
    
    
    return(list("df"=df, "scale"=scale, "acc_rate"=acc_rate))
}