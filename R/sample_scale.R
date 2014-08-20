sample_scale <- function(c, m, z, U_xi, U_psi, 
                         U_Sigma, U_df, ltn, weights, scale){
    
    
    
#     prior_df <- function(x, log=FALSE){
#         if(log){
#             y <- log(x/100)-x/10
#         }else{
#             y <- x/100*exp(-x/10)
#         }
#         return(y)
#     }
    
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
    prior_df <- function(x, log=FALSE){
        if(log){
            y <- 1/2*(log(x)-log(x+3) +log(trigamma(x/2)-trigamma((x+1)/2)-2*(x+3)/(x*(x+1)^2)))
        }else{
            y <- sqrt(x/(x+3)*(trigamma(x/2)-trigamma((x+1)/2)-2*(x+3)/(x*(x+1)^2)))
        }
        return(y)
    }
    
    df_new <- list()
    c_df <- 4
    fullCl_ind <- which(m!=0)
    fullCl <- length(fullCl_ind)
    
    
    U_xi_list <- lapply(fullCl_ind, function(j) U_xi[, j])
    U_psi_list <- lapply(fullCl_ind, function(j) U_psi[, j])
    U_Sigma_list <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
    U_df_list <- lapply(fullCl_ind, function(j) U_df[j])
    
    df_new <- list()
    for(j in 1:fullCl){
        df <- U_df_list[[j]]
        df_new[[j]] <- 1+exp(runif(1, min=log(df-1)-c_df, max = log(df-1)+c_df))
    }
    
    
    likelihood_old <- apply(X=mvstpdf(z, xi=U_xi_list, 
                                      sigma=U_Sigma_list, 
                                      psi=U_psi_list, df=U_df_list), 
                            MARGIN=1, FUN="*",y=weights[fullCl_ind])
    likelihood_new <- apply(X=mvstpdf(z, xi=U_xi_list, 
                                      sigma=U_Sigma_list, 
                                      psi=U_psi_list, df=df_new), 
                            MARGIN=1, FUN="*",y=weights[fullCl_ind])
    
    u <- runif(fullCl)
    if(fullCl>1){
        for(j in 1:fullCl){
            prob_new <- exp(sum(log(colSums(rbind(likelihood_old[-j,],likelihood_new[j,]))) + prior_df(df_new[[j]], log=TRUE) + log(df_new[[j]]-1))
                            -sum(log(colSums(likelihood_old)) + prior_df(U_df[[j]], log=TRUE) + log(U_df[[j]]-1)))
            if(prob_new==Inf){
                prob_new <- 1
            }else{ 
                prob_new <- min(1, prob_new)
            }
            if (u[j]<prob_new){
                U_df_list[[j]] <- df_new[[j]]
            }
            
            obs_j <- which(c==fullCl_ind[j])
            
            eps <- z[,obs_j] - U_xi_list[[j]] - sapply(X=ltn[obs_j], FUN=function(x){x*U_psi_list[[j]]})
            tra <- apply(X=eps, MARGIN=2, FUN=function(v){sum(diag(tcrossprod(v)%*%solve(U_Sigma_list[[j]])))})
            
            scale[obs_j] <- rgamma(length(obs_j), 
                                   shape=(U_df_list[[j]] + nrow(z) + 1)/2, 
                                   rate=(U_df_list[[j]] + ltn[obs_j]^2 + tra)/2)
        }
    }else{
        j <- 1
        prob_new <- exp(sum(log(likelihood_new) + prior_df(df_new[[j]], log=TRUE) + log(df_new[[j]]-1))
                        -sum(log(likelihood_old) + prior_df(U_df[[j]], log=TRUE) +log(U_df[[j]]-1)))
        if(prob_new==Inf){
            prob_new <- 1
        }else{ 
            prob_new <- min(1, prob_new)
        }
        if (u[j]<prob_new){
            U_df_list[[j]] <- df_new[[j]]
        }
        
        obs_j <- which(c==fullCl_ind[j])
        
        eps <- z[,obs_j] - U_xi_list[[j]] - sapply(X=ltn[obs_j], FUN=function(x){x*U_psi_list[[j]]})
        tra <- apply(X=eps, MARGIN=2, FUN=function(v){sum(diag(tcrossprod(v)%*%solve(U_Sigma_list[[j]])))})
        
        scale[obs_j] <- rgamma(length(obs_j), 
                               shape=(U_df_list[[j]] + nrow(z) + 1)/2, 
                               rate=(U_df_list[[j]] + ltn[obs_j]^2 + tra)/2)
    }
    
    
    return(list("df"=U_df_list, "scale"=scale))
}