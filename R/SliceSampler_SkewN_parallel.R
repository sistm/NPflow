sliceSampler_SkewN_parallel <- function(c, m, alpha, z, hyperG0, U_xi, U_psi, U_Sigma){
    
    maxCl <- length(m) #maximum number of clusters
    ind <- which(m!=0) # indexes of non empty clusters
    
    # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
    # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1)    # and sample the rest of the weigth for potential new clusters:
    # temp_{K+1} ~ Gamma(alpha, 1)
    # then renormalise temp
    w <- numeric(maxCl)
    temp <- rgamma(n=(length(ind)+1), shape=c(m[ind], alpha), scale = 1)
    temp_norm <- temp/sum(temp)
    w[ind] <- temp_norm[-length(temp_norm)]
    R <- temp_norm[length(temp_norm)] 
    #R is the rest, i.e. the weight for potential new clusters
    
    # Sample the latent u
    u  <- runif(maxCl)*w[c]
    u_star <- min(u)
    
    # Sample the remaining weights that are needed with stick-breaking
    # i.e. the new clusters
    ind_new <- which(m==0) # potential new clusters
    if(length(ind_new)>0){
        t <- 0 # the number of new non empty clusters
        while(R>u_star && (t<length(ind_new))){ 
            # sum(w)<1-min(u) <=> R>min(u) car R=1-sum(w)
            t <- t+1
            beta_temp <- rbeta(n=1, shape1=1, shape2=alpha)
            # weight of the new cluster
            w[ind_new[t]] <- R*beta_temp
            R <- R * (1-beta_temp) # remaining weight
        }
        ind_new <- ind_new[1:t]
        
        # Sample the centers and spread of each new cluster from prior
        for (i in 1:t){
            NNiW <- nniw_rnd(hyperG0)
            #TODO
            U_xi[, ind_new[i]] <- NNiW[["xi"]]
            U_psi[, ind_new[i]] <- NNiW[["psi"]]
            U_Sigma[, , ind_new[i]] <- NNiW[["S"]]
        }
    }
    
    fullCl_ind <- which(w != 0)
    nb_fullCl_ind <- length(fullCl_ind)
    # calcul de la vraisemblance pour chaque données pour chaque clusters
    # assignation de chaque données à 1 cluster
     # likelihood of belonging to each cluster 
    
    
    browser()
    #TODO
#     U_xi_list <- lapply(fullCl_ind, function(j) U_xi[, j])
#     U_psi_list <- lapply(fullCl_ind, function(j) U_psi[, j])
#     U_Sigma_list <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
#     l <- apply(X=mvsnpdf(z, xi=U_xi_list, sigma=U_Sigma_list, psi=U_psi_list), MARGIN=1, FUN="*",y=w[fullCl_ind])          
#     c <- apply(X=l, MARGIN=2, FUN=which.max)
    
    c <- foreach(i=1:maxCl, .combine='c')%dopar%{
        l <- numeric(nb_fullCl_ind)
        for(j in fullCl_ind){
            l[j] <- mvsnpdf(x = matrix(z[,i], ncol= 1, nrow=length(z[,i])) , 
                            xi = U_xi[, j], 
                            sigma = U_Sigma[, , j],
                            psi = U_psi[,j]
            
            )*w[j]            
        }
        c_par <- which.max(l)
    }
    
    m_new <- numeric(maxCl) # number of observations in each cluster
    m_new[unique(c)] <- table(c)
    
    ltn <- numeric(maxCl) # latent truncated normal variables
    for (k in which(m_new!=0)){
        obs_k <- which(c==k)
        siginv <- solve(U_Sigma[, , k])
        psi <- U_psi[,k]
        A_k <-  1/(1 + (crossprod(psi, siginv)%*%psi))
        a_ik <- (tcrossprod(A_k, psi)%*%siginv%*%(z[,obs_k]-U_xi[,k]))
        ltn[obs_k] <- rtruncnorm(length(obs_k), a=0, b=Inf, mean = a_ik, sd = sqrt(A_k))
    }
    
    return(list("c"=c, "m"=m_new, "weights"=w, "latentTrunc"=ltn))
}