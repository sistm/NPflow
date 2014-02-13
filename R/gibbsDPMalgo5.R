#'Gibbs Sapling DPM Alghorithm 4: slice sampling
#'with a prior on alpha
#'
#'@param z data
#'
#'@param hyperG0 prior mixing distribution
#'
#'@param a shape parameter f the Gamma hyperprior
#'
#'@param b scale parameter f the Gamma hyperprior
#'
#'@param N number of MCMC iterations
#'
#'@author Boris Hejblum
#'
#'@export gibbsDPMalgo5
#'
#'@examples
#' rm(list=ls())
#' #Number of data
#' n <- 2000
#' set.seed(1234)
#' #set.seed(4321)
#' 
#' # Sample data
#' m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2))
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
#'  cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'  
#'  # Set parameters of G0
#'  hyperG0 <- list()
#'  hyperG0[["mu"]] <- c(0,0)
#'  hyperG0[["kappa"]] <- 1
#'  hyperG0[["nu"]] <- 4
#'  hyperG0[["lambda"]] <- diag(2)
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.1
#'  b <- 0.001
#'  plot(density(rgamma(n=5000, a, 1/b)))
#'  # Number of iterations
#'  N <- 15
#'  
#'  # do some plots
#'  doPlot <- TRUE 
#'  nbclust_init <- 30
#'  
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  MCMCsample <- gibbsDPMalgo5(z, hyperG0, a, b, N, doPlot, nbclust_init)
#'  
#'  alpha_m <- mean(MCMCsample$alpha[floor(length(MCMCsample$alpha)/2):length(MCMCsample$alpha)])
#'  sum(alpha_m/(alpha_m+1:n-1)) # expected number of clusters
#'  
#'  plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster)
#'
#'
gibbsDPMalgo5 <- function (z, hyperG0, a, b, N, doPlot=TRUE, nbclust_init=30){
    
    if(doPlot){library(ggplot2)}
    
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    
    # U_SS is a list where each U_SS[[k]] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    #store U_SS :
    U_SS_list <- list()
    #store clustering :
    c_list <- list()
    
    m <- numeric(n) # number of obs in each clusters
    c <-numeric(n)
     # initial number of clusters
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    if(ncol(z)<nbclust_init){       
        for (k in 1:n){
            c[k] <- k
            U_SS[[c[k]]] <- update_SS(z=z[, k], S=hyperG0)
            NiW <- normalinvwishrnd(U_SS[[c[k]]])
            U_mu[, c[k]] <- NiW[["mu"]]
            U_Sigma[, , c[k]] <- NiW[["S"]]
            m[c[k]] <- m[c[k]]+1
        }
    } else{
        c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
        for (k in unique(c)){
            obs_k <- which(c==k)
            U_SS[[k]] <- update_SS(z=z[, obs_k], S=hyperG0)
            NiW <- normalinvwishrnd(U_SS[[k]])
            U_mu[, k] <- NiW[["mu"]]
            U_Sigma[, , k] <- NiW[["S"]]
            m[k] <- m[k]+1
        }
    }
    
    
    
    nbClust <- length(unique(c))
    alpha <- log(n)
    
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    
    cat(i, "/", N, " samplings\n", sep="")
    if(doPlot){
        plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, 
                 m=m, c=c, i=i, alpha=alpha[length(alpha)])
    }
    
    
    for(i in 2:N){
        
        alpha <- c(alpha,
                   sample_alpha(alpha_old=alpha[length(alpha)], n=n, 
                                K=nbClust, a=a, b=b)
        )
        slice <- slice_sample(c=c, m=m, alpha=alpha[length(alpha)], 
                              z=z, hyperG0=hyperG0, 
                              U_mu=U_mu, U_Sigma=U_Sigma)
        m <- slice[["m"]]
        c <- slice[["c"]]
        U_mu <- slice[["U_mu"]]
        U_Sigma <- slice[["U_Sigma"]]
        
        
        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in fullCl){
            obs_j <- which(c==j)
            if(j > length(U_SS) || is.null(U_SS[[j]])){
                U_SS[[j]] <- update_SS(z=z[, obs_j], S=hyperG0)
            } else{
                U_SS[[j]] <- update_SS(z[,obs_j], S=U_SS[[j]])
            }
            NiW <- normalinvwishrnd(U_SS[[j]])
            U_mu[, j] <- NiW[["mu"]]
            U_Sigma[, , j] <- NiW[["S"]]
        }
        
        
        U_SS_list[[i]] <- U_SS[which(m!=0)]
        c_list[[i]] <- c
        
        cat(i, "/", N, " samplings\n", sep="")
        if(doPlot){
            plot_DPM(z, U_mu, U_Sigma, m, c, i, alpha=alpha[length(alpha)])
        }
    
    }
    
    return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
                "partition" = m, "alpha"=alpha, "U_SS_list"=U_SS_list,
                "c_list" = c_list))
}






