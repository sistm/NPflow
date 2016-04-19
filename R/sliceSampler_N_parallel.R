#'@keywords internal
#'@author Boris Hejblum
#'@importFrom stats rbeta rgamma runif
sliceSampler_N_parallel <- function(Ncpus, c, m, alpha, z, hyperG0, U_mu, U_Sigma, diagVar, parallel_index){

    maxCl <- length(m) #maximum number of clusters
    ind <- which(m!=0) #indexes of non empty clusters
    nb_fullCl <- length(ind) #number of non empty clusters

    # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
    # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1)    # and sample the rest of the weigth for potential new clusters:
    # temp_{K+1} ~ Gamma(alpha, 1)
    # then renormalise temp
    w <- numeric(maxCl)
    temp <- stats::rgamma(n=(nb_fullCl+1), shape=c(m[ind], alpha), scale = 1)
    temp_norm <- temp/sum(temp)
    w[ind] <- temp_norm[-(nb_fullCl+1)]
    R <- temp_norm[(nb_fullCl+1)]
    #R is the rest, i.e. the weight for potential new clusters


    # Sample the latent u
    u  <- stats::runif(maxCl)*w[c]
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
            beta_temp <- stats::rbeta(n=1, shape1=1, shape2=alpha)
            # weight of the new cluster
            w[ind_new[t]] <- R*beta_temp
            R <- R * (1-beta_temp) # remaining weight
        }

        if(t>0){
            ind_new <- ind_new[1:t]
            # Sample the centers and spread of each new cluster from prior
            for (i in 1:t){
                NiW <- rNiW(hyperG0, diagVar)
                U_mu[[as.character(ind_new[i])]] <- NiW[["mu"]]
                U_Sigma[[as.character(ind_new[i])]] <- NiW[["S"]]
            }
        }
    }

    fullCl_ind <- which(w != 0)

    # likelihood of belonging to each cluster computation
    # sampling clusters
    if(length(fullCl_ind)>1){
        U_mu_full <- sapply(fullCl_ind, function(j) U_mu[, j])
        U_Sigma_full <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])

        c <- foreach::"%dopar%"(foreach::foreach(i=1:Ncpus, .combine='c'),
                                {
            l <- mmvnpdfC(x=z[, parallel_index[[i]]], mean=U_mu_full, varcovM=U_Sigma_full, Log = FALSE)
            u_mat <- t(sapply(w[fullCl_ind], function(x){as.numeric(u[parallel_index[[i]]] < x)}))
            prob_mat <- u_mat * l

            #fast C++ code
            c <- fullCl_ind[sampleClassC(prob_mat)]
            #         #slow C++ code
            #         c <- fullCl_ind[sampleClassC_bis(prob_mat)]
            #         #vectorized R code
            #         c <- fullCl_ind[apply(X= prob_mat, MARGIN=2, FUN=function(v){match(1,rmultinom(n=1, size=1, prob=v))})]
            #         #alternative implementation:
            #         prob_colsum <- colSums(prob_mat)
            #         prob_norm <- apply(X=prob_mat, MARGIN=1, FUN=function(r){r/prob_colsum})
            #         c <- fullCl_ind[apply(X=prob_norm, MARGIN=1, FUN=function(r){match(TRUE,stats::runif(1) <cumsum(r))})]
        })
    }else{
        c <- rep(fullCl_ind, maxCl)
    }

    m_new <- numeric(maxCl) # number of observations in each cluster
    m_new[unique(c)] <- table(c)

    return(list("c"=c, "m"=m_new, "weights"=w))
}