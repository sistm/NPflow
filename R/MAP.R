#'Maximum a Posteriori from MCMC samples of Dirichlet Process Mixture of Skew Normal
#'
#'Hyperprior on the concentration parameter of the DP
#'
#'@param s an DPM object 
#'
#'@param from an integer giving the first iteration from which MAP is 
#'to be looked for
#'
#'
#'
#'
MAP <- function(s, from=1){
    
    logposts <- unlist(lapply(s$logposterior_list, sum))
    N_iter <- length(logposts)
    map_ind <- from - 1 + which.max(logposts[c(from:N_iter)])
    map <- MCMCsample_sn$c_list[[map_ind]]
    
    
    clusters <- s$c_list[[map_ind]]
    freq <- table(clusters)
    
    U_SS <- s$U_SS_list[[map_ind]]
    names(U_SS) <- names(freq)
    if(!is.null(U_SS[[1]][["xi"]])){
        U_xi <- lapply(U_SS, "[[", "xi")
        U_psi <- lapply(U_SS, "[[", "psi")
    }
    U_Sigma <- lapply(U_SS, "[[", "S")
    
    partition <- numeric(length(s$partition))
    partition[as.integer(names(freq))] <- freq
    
    return(list("iteration_MAP" = map_ind,
                "nb_iterations" = N_iter,
                "clusters" = s$c_list[[map_ind]], 
                "U_SS" = U_SS,
                "U_xi"=U_xi,
                "U_psi"=U_psi,
                "U_Sigma"=U_Sigma,
                "partition"=partition,
                "alpha"=s$alpha[map_ind],
                "logposterior"=logposts[map_ind]))
    
}
    