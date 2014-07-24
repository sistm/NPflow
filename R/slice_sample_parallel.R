slice_sample_parallel <- function(c, m, alpha, z, hyperG0, U_mu, U_Sigma){
    
    maxCl <- length(m) #maximum number of clusters
    ind <- which(m!=0) #indexes of non empty clusters
    nb_fullCl <- length(ind) #number of non empty clusters
    
    # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
    # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1)    # and sample the rest of the weigth for potential new clusters:
    # temp_{K+1} ~ Gamma(alpha, 1)
    # then renormalise temp
    w <- numeric(maxCl)
    temp <- rgamma(n=(nb_fullCl+1), shape=c(m[ind], alpha), scale = 1)
    temp_norm <- temp/sum(temp)
    w[ind] <- temp_norm[-(nb_fullCl+1)]
    R <- temp_norm[(nb_fullCl+1)] 
    #R is the rest, i.e. the weight for potential new clusters
    
    
    # Sample the latent u
    u  <- runif(maxCl)*w[c]
    u_star <- min(u)
    
    # Sample the remaining weights that are needed with stick-breaking
    # i.e. the new clusters
    if(nb_fullCl < maxCl){
        ind_new <- c((nb_fullCl+1):maxCl) # potential new clusters
    }else{
        ind_new <- integer(0)
    }
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
        
        if(t>0){
            ind_new <- ind_new[1:t]
            # Sample the centers and spread of each new cluster from prior
            for (i in 1:t){
                NiW <- rNiW(hyperG0)
                U_mu[[as.character(ind_new[i])]] <- NiW[["mu"]]
                U_Sigma[[as.character(ind_new[i])]] <- NiW[["S"]]
            }
        }
    }
    nb_fullCl <- nb_fullCl + t
    fullCl_ind <- which(w != 0)
    # calcul de la vraisemblance pour chaque données pour chaque clusters
    # assignation de chaque données à 1 cluster
    # likelihood of belonging to each cluster 
    c <- foreach(i=1:maxCl, .combine='c')%dopar%{
        l <- numeric(nb_fullCl)
        for (j in fullCl_ind){
            k <- as.character(j)
            l[j] <- mvnpdf(x = matrix(z[,i], ncol= 1, nrow=length(z[,i])) , 
                           mean = U_mu[[k]], 
                           varcovM = U_Sigma[[k]])*w[j]  
        }
        c_par <- which.max(l)
    }
    
    m_new <- numeric(maxCl) # number of observations in each cluster
    m_new[unique(c)] <- table(c)
    
    return(list("c"=c, "m"=m_new, "weights"=w))
}