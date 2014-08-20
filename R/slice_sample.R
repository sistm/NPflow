slice_sample <- function(c, m, alpha, z, hyperG0, U_mu, U_Sigma){
    
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
            NiW <- rNiW(hyperG0)
            U_mu[, ind_new[i]] <- NiW[["mu"]]
            U_Sigma[, , ind_new[i]] <- NiW[["S"]]
        }
    }
    
    fullCl_ind <- which(w != 0)
    # calcul de la vraisemblance pour chaque données pour chaque clusters
    # assignation de chaque données à 1 cluster
    if(length(fullCl_ind)>1){
        U_mu_list <- lapply(fullCl_ind, function(j) U_mu[, j])
        U_Sigma_list <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
        l <- mvnpdf(z, mean=U_mu_list, varcovM=U_Sigma_list)
        u_mat <- apply(X=t(sapply(u, function(x){x < w[fullCl_ind]})), MARGIN=2, FUN= as.numeric)
        prob_mat <- u_mat * l
        c <- apply(X= prob_mat, MARGIN=1, FUN=function(v){which(rmultinom(n=1, size=1, prob=v)==1)})
    }else{
        c <- rep(fullCl_ind, maxCl)
    }
    m_new <- numeric(maxCl) # number of observations in each cluster
    m_new[unique(c)] <- table(c)[as.character(unique(c))]
    
    # non vectorized code for cluster allocation:
    #     l <- numeric(length(fullCl_ind)) # likelihood of belonging to each cluster
    #     m_new <- numeric(maxCl) # number of observations in each cluster
    #     for(i in 1:maxCl){
    #         for (j in fullCl_ind){
    #             l[j] <- mvnpdf(x = matrix(z[,i], ncol= 1, nrow=length(z[,i])) , 
    #                            mean = U_mu[, j], 
    #                            varcovM = U_Sigma[, , j])*w[j]  
    #         }
    #         c[i] <- which.max(l)
    #         m_new[c[i]] <- m_new[c[i]] + 1
    #     }
    
    return(list("c"=c, "m"=m_new, "weights"=w))
}