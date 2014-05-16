#'Slice Sampling of Dirichlet Process Mixture of skew Normals
#'
#'Hyperprior on the concentration parameter of the DP
#'
#'@param z data matrix \code{d x n} with \code{d} dimensions in rows 
#'and \code{n} observations in columns.
#'
#'@param hyperG0 prior mixing distribution.
#'
#'@param a shape hyperparameter of the Gamma prior 
#'on the parameter of the Dirichlet Process.
#'
#'@param b scale hyperparameter of the Gamma prior 
#'on the parameter of the Dirichlet Process.
#'
#'@param N number of MCMC iterations.
#'
#'@param doPlot logical flag indicating wether to plot MCMC iteration or not.
#'Default to \code{TRUE}.
#'
#'@param nbclust_init number of clusters at initialisation. 
#'Default to 30 (or less if there are less than 30 observations).
#'
#'@author Boris Hejblum
#'
#'@export DPM_GibbsSampler_SkewN_test1K
#'
#'@examples
#' rm(list=ls())
#' library(ggplot2)
#' #Number of data
#' n <- 1000
#' #n <- 2000
#' set.seed(1234)
#' #set.seed(123)
#' #set.seed(4321)
#' 

#' 
#'# Sample data
#' d <- 2
#' xi <- matrix(nrow=d, ncol=1, c(-1.5, 1))
#' psi <- matrix(nrow=d, ncol=1, c(3, -2))
#' p <- 1 # frequence des clusters
#' 
#' sdev <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))


#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + matrix(rnorm(d, mean = 0, sd = 1), nrow=1, ncol=d)%*%sdev
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'  

#'  # Set parameters of G0
#'  hyperG0 <- list()
#'  hyperG0[["b_xi"]] <- rep(0,d) # c(-1.5, 1) for an informative prior instead...
#'  hyperG0[["b_psi"]] <- rep(0,d) # c(3, -2) for an informative prior instead...
#'  hyperG0[["kappa"]] <- 0.001
#'  hyperG0[["D_xi"]] <- 1000 # weak prior
#'  hyperG0[["D_psi"]] <- 1000 #0.00001 for a strong prior
#'  hyperG0[["nu"]] <- 0.1
#'  hyperG0[["lambda"]] <- diag(d)/10 # diag(d)*0.3 for an informative prior instead...
#'  
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'  
#'  # do some plots
#'  doPlot <- TRUE 
#'  nbclust_init <- 30
#'  
#'  
#'  
#'  ## Data
#'  ########
#'  p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y)) 
#'        + geom_point()
#'        + ggtitle("Simple example in 2d data")
#'        +xlab("D1")
#'        +ylab("D2")
#'        +theme_bw())
#'  p
#'  
#'  
#'  ## alpha priors plots
#'  #####################
#'  prioralpha <- data.frame("alpha"=rgamma(n=5000, shape=a, scale=1/b), 
#'                          "distribution" =factor(rep("prior",5000), 
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(prioralpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="red")
#'        + ggtitle(paste("Prior distribution on alpha: Gamma(", a, 
#'                  ",", b, ")\n", sep=""))
#'       )
#'  p
#'  
#'  
#'  
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'  MCMCsample_sn <- DPM_GibbsSampler_SkewN_test1K(z, hyperG0, a, b, N=200, doPlot, nbclust_init, plotevery=10)
#'  MCMCsample_sn$U_xi
#'  MCMCsample_sn$U_psi
#'
#'
#'
#'

DPM_GibbsSampler_SkewN_test1K <- function (z, hyperG0, a, b, N, doPlot=TRUE, 
                                           nbclust_init=30, plotevery=1, ...){
    
    if(doPlot){library(ggplot2)}
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_xi <- numeric(p)
    U_psi <- numeric(p)
    U_Sigma = matrix(0, ncol=p, nrow=p)
    
    # U_SS is a list where each U_SS[[k]] contains the sufficient
    # statistics associated to cluster k
    
    #store U_SS :
    U_SS_list <- list()
    #store clustering :
    c_list <- list()
    #store sliced weights
    weights_list <- list()
    
    #store log posterior probability
    logposterior_list <- list()
    
    m <- numeric(n) # number of obs in each clusters
    c <- numeric(n) # cluster label of ech observation
    ltn <- rtruncnorm(n, a=0, b=Inf, mean=0, sd=1) # latent truncated normal
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1      
    U_SS <- update_SSsn(z, S=hyperG0, ltn=ltn)
    NNiW <- nniw_rnd(U_SS)
    U_xi <- NNiW[["xi"]]
    U_psi <- NNiW[["psi"]]
    U_Sigma <- NNiW[["S"]]
    
    alpha <- c(log(n))
    
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[1]] <- numeric(length(m))
    weights_list[[1]][unique(c)] <- table(c)/length(c)
    
    logposterior_list[[i]] <- logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_SS[["B"]], 
                                                 hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
    
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    
    if(doPlot){
        plot_DPMsn_test1K(z=z, U_xi=U_xi, U_psi=U_psi, U_Sigma=U_Sigma, 
                          m=m, c=c, i=i, alpha=alpha[i], U_SS=U_SS, ellipses=F, ...)
    }else{
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    if(N>1){
        for(i in 2:N){
            nbClust <- 1
            
            alpha <- c(alpha,
                       sample_alpha(alpha_old=alpha[i-1], n=n, 
                                    K=nbClust, a=a, b=b)
            )
            slice <- sliceSampler_SkewN_test1K(c=c, m=m, alpha=alpha[i], 
                                               z=z, hyperG0=hyperG0, 
                                               U_xi=U_xi, U_psi=U_psi, U_Sigma=U_Sigma)
            
            ltn <- slice[["latentTrunc"]]
            
            # Update locations
            
            U_SS <- update_SSsn(z=z, S=hyperG0,  ltn=ltn)
            NNiW <- nniw_rnd(U_SS)
            U_xi <- NNiW[["xi"]]
            U_psi <- NNiW[["psi"]]
            U_Sigma <- NNiW[["S"]]
            
            
            U_SS_list[[i]] <- U_SS
            c_list[[i]] <- c
            
            logposterior_list[[i]] <- logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_SS[["B"]], 
                                                         hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
            
            cat(i, "/", N, " samplings:\n", sep="")
            cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
            
            if(doPlot && i/plotevery==floor(i/plotevery)){
                plot_DPMsn_test1K(z=z, U_xi=U_xi, U_psi=U_psi, U_Sigma=U_Sigma, m=m, c=c, i=i,
                                  alpha=alpha[i], U_SS=U_SS, ellipses=F, ...)
            }else{
                cl2print <- unique(c)
                cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
            }
            
        }
    }
    
    return(list("clusters" = c, "U_xi" = U_xi, "U_psi" = U_psi, "U_Sigma" = U_Sigma, 
                "partition" = m, "alpha"=alpha, "U_SS_list"=U_SS_list,
                "c_list" = c_list, "weights_list"=weights_list, 
                "logposterior_list"=logposterior_list))
}






