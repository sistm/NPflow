#'Gibbs Sapling DPM Alghorithm 4: slice sampling
#'
#'@param z
#'
#'@param hyperG0
#'
#'@param alpha
#'
#'@param N
#'
#'@author Francois Caron
#'
#'@export gibbsDPMalgo4
#'
#'@examples
#' rm(list=ls())
#' #Number of data
#' n <- 100
#' #set.seed(1231)
#' 
#' # Sample data
#' m <- matrix(nrow=2, ncol=4, c(-1, 1, 0, 2, 1, -2, -1, -2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' 
#' library(expm)
#' s <- array(dim=c(2,2,4))
#' s[, ,1] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.1))
#' s[, ,2] <- matrix(nrow=2, ncol=2, c(0.01, 0, 0, 0.1))
#' s[, ,3] <- matrix(nrow=2, ncol=2, c(0.1, 0.08, 0.08, 0.1))
#' s[, ,4] <- .1*diag(2)
#' c <- rep(0,n)
#' z <- matrix(0, nrow=2, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- m[, c[k]] + sqrtm(s[, , c[k]])%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
#' }
#'  
#'  # Set parameters of G0
#'  hyperG0 <- list()
#'  hyperG0[["mu"]] <- c(0,0)
#'  hyperG0[["kappa"]] <- 1
#'  hyperG0[["nu"]] <- 4
#'  hyperG0[["lambda"]] <- diag(2)
#'  # Scale parameter of DPM
#'  alpha <- 3
#'  # Number of iterations
#'  N <- 10 
#'  # do some plots
#'  doPlot <- TRUE 
#'  
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  GibSample <- gibbsDPMalgo4(z, hyperG0, alpha, N, doPlot)
#'
#'
gibbsDPMalgo4 <- function (z, hyperG0, alpha, N, doPlot=TRUE){
    
    if(doPlot){library(ggplot2)}
    
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    
    # U_SS is a list where each U_SS[k] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    m <- numeric(n) # number of obs in each clusters
    c <-numeric(n)
    
    # Initialisation: each observation is assigned to a different cluster----
    
    i <- 1
    for (k in 1:n){
        c[k] <- k
        U_SS[[c[k]]] <- update_SS(z=z[, k], S=hyperG0)
        NiW <- normalinvwishrnd(U_SS[[c[k]]])
        U_mu[, c[k]] <- NiW[["mu"]]
        U_Sigma[, , c[k]] <- NiW[["S"]]
        m[c[k]] <- m[c[k]]+1
    }
    
    cat(i, "/", N, " samplings\n", sep="")
    if(doPlot){
        plot_DPM4(z, U_mu, m, c, i)
    }
    
    
    for(i in 2:N){
        
        slice <- slice_sample(c=c, m=m, alpha=alpha, z=z, 
                              hyperG0=hyperG0, 
                              U_mu=U_mu, U_Sigma=U_Sigma)
        m <- slice[["m"]]
        c <- slice[["c"]]
        U_mu <- slice[["U_mu"]]
        U_Sigma <- slice[["U_Sigma"]]
        
        # Update cluster locations U
        ind <- which(m!=0)
        
        for(j in 1:length(ind)){
            NiW <- normalinvwishrnd(U_SS[[ind[j]]])
            U_mu[, ind[j]] <- NiW[["mu"]]
            U_Sigma[, , ind[j]] <- NiW[["S"]]
        }
        
        cat(i, "/", N, " samplings\n", sep="")
        if(doPlot){
            plot_DPM4(z, U_mu, m, c, i)
        }
    }
    return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
                "partition" = m))
}








# Subfunctions ----

slice_sample <- function(c, m, alpha, z, hyperG0, U_mu, U_Sigma){
    
    maxCl <- length(m) #maximum number of clusters
    ind <- unique(c) # non empty clusters
    fullCl <- which(m!=0) # indexes of non empty clusters
    r <- sum(m)
    
    # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
    # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1), temp_{K+1} ~ Gamma(gamma, 1)
    #renormalisation of temp
    w <- numeric(maxCl)
    temp <- rgamma(n=(length(ind)+1), shape=c(m[ind], alpha), scale = 1)
    #temp = gamrnd([m(ind); gamma], 1);
    temp_norm <- temp/sum(temp)
    w[ind] <- temp_norm[-length(temp_norm)]
    R <- temp_norm[length(temp_norm)] #the rest of the wights
    
    
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
            NiW <- normalinvwishrnd(hyperG0)
            U_mu[, ind_new[i]] <- NiW[["mu"]]
            U_Sigma[, , ind_new[i]] <- NiW[["S"]]
        }
    }
    
    # calcul de la vraisemblance pour chaque données pour chaque clusters
    # assignation de chaque données à 1 cluster
    l <- numeric(length(fullCl)) # likelihood of belonging to each cluster 
    m <- numeric(maxCl) # number of observations in each cluster
    for(i in 1:maxCl){
        for (j in 1:length(fullCl)){
            l[j] <- mvnpdf(x = matrix(z[,i], nrow= 1, ncol=length(z[,i])) , 
                             mean = U_mu[, fullCl[j]], 
                             varcovM = U_Sigma[, , fullCl[j]])*w[fullCl[j]]  
        }
        c[i] <- which.max(l)
        m[c[i]] <- m[c[i]] + 1
    }
    
    return(list("c"=c, "m"=m, "U_mu"=U_mu, "U_Sigma"=U_Sigma))
}





plot_DPM4 <- function(z, U_mu, m, c, i){
    fullCl <- which(m!=0)
    U_mu2plot <- U_mu[, fullCl]    
    zClusters <- as.factor(c)
    levels(zClusters) <- as.character(1:length(levels(zClusters)))
    z2plot <- cbind.data.frame("X"=z[1,],"Y"=z[2,],"Cluster"=zClusters)
    U2plot <- cbind.data.frame("X"=U_mu2plot[1,],"Y"=U_mu2plot[2,],"Cluster"=factor(1:dim(U_mu2plot)[2]))
    p <- (ggplot(z2plot) 
          + geom_point(aes(x=X, y=Y, col=Cluster), data=z2plot) 
          + geom_point(aes(x=X, y=Y, col=Cluster), data=U2plot, shape="X", size=5)
          + ggtitle(paste("Gibbs sampling for DPM - algo 2\nIteration", i))
    )
    print(p)
}
