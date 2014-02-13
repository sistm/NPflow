#'Gibbs Sapling DPM Alghorithm 4: slice sampling
#'
#'@param z data
#'
#'@param hyperG0 prior mixing distribution
#'
#'@param alpha parameter of the Dirichlet Process
#'
#'@param N number of MCMC iterations
#'
#'@author Boris Hejblum
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
#'  alpha <- 4
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
    ninit_clust <- 50
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    if(ncol(z)<ninit_clust){       
        for (k in 1:n){
            c[k] <- k
            U_SS[[c[k]]] <- update_SS(z=z[, k], S=hyperG0)
            NiW <- normalinvwishrnd(U_SS[[c[k]]])
            U_mu[, c[k]] <- NiW[["mu"]]
            U_Sigma[, , c[k]] <- NiW[["S"]]
            m[c[k]] <- m[c[k]]+1
        }
    } else{
        c <- sample(x=1:ninit_clust, size=n, replace=TRUE)
        for (k in unique(c)){
            obs_k <- which(c==k)
            U_SS[[k]] <- update_SS(z=z[, obs_k], S=hyperG0)
            NiW <- normalinvwishrnd(U_SS[[k]])
            U_mu[, k] <- NiW[["mu"]]
            U_Sigma[, , k] <- NiW[["S"]]
            m[k] <- m[k]+1
        }
    }
        
    cat(i, "/", N, " samplings\n", sep="")
    if(doPlot){
        plot_DPM(z, U_mu, m, c, i)
    }
    
    
    for(i in 2:N){
        
        slice <- slice_sample(c=c, m=m, alpha=alpha, z=z, 
                              hyperG0=hyperG0, 
                              U_mu=U_mu, U_Sigma=U_Sigma)
        m <- slice[["m"]]
        c <- slice[["c"]]
        U_mu <- slice[["U_mu"]]
        U_Sigma <- slice[["U_Sigma"]]
        
        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in 1:length(fullCl)){
            obs_j <- which(c==fullCl[j])
            if(is.null(U_SS[[fullCl[j]]])){
                U_SS[[fullCl[j]]] <- update_SS(z=z[, obs_j], S=hyperG0)
            } else{
                U_SS[[fullCl[j]]] <- update_SS(z[,obs_j], 
                                               S=U_SS[[fullCl[j]]])
            }
            NiW <- normalinvwishrnd(U_SS[[fullCl[j]]])
            U_mu[, fullCl[j]] <- NiW[["mu"]]
            U_Sigma[, , fullCl[j]] <- NiW[["S"]]
        }
        
        cat(i, "/", N, " samplings\n", sep="")
        if(doPlot){
            plot_DPM(z, U_mu, m, c, i)
        }
    }
    return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
                "partition" = m))
}




