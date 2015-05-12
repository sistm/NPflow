#'Slice Sampling of Dirichlet Process Mixture of skew Normals
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
#'@param diagVar logical flag indicating wether the variance of each cluster is 
#'estimated as a diagonal matrix, or as a full matrix. 
#'Default is \code{TRUE} (diagonal variance).
#'
#'@param nbclust_init number of clusters at initialisation. 
#'Default to 30 (or less if there are less than 30 observations).
#'
#'@return a object of class \code{DPMclust} with the following attributes: 
#'  \itemize{
#'      \item{\code{mcmc_partitions}:}{a list of length \code{N}. Each
#'       element \code{mcmc_partitions[n]} is a vector of length 
#'       \code{n} giving the partition of the \code{n} observations.}
#'      \item{\code{alpha}:}{ a vector of length \code{N}. \code{cost[j]} is the cost 
#' associated to partition \code{c[[j]]}}
#'      \item{\code{weights_list}:}{}
#'      \item{\code{logposterior_list}:}{}
#'      \item{\code{data}:}{the data matrix \code{d x n} with \code{d} dimensions in rows 
#'and \code{n} observations in columns.}
#'      \item{\code{nb_mcmcit}:}{the number of MCMC itertations}
#'  }
#'
#'@author Boris Hejblum
#'
#'@export
#'
#'@examples
#' rm(list=ls())
#' library(ggplot2)
#' 
#' #Number of data
#' n <- 2000
#' set.seed(123)
#' set.seed(1234)
#' #set.seed(4321)
#' 
#' 
#' d <- 2
#' ncl <- 4
#' 
#' # Sample data
#' 
#' sdev <- array(dim=c(d,d,ncl))
#' 
#' xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1, 1.5, 1, 1.5, -2, -2, -2))
#' psi <- matrix(nrow=d, ncol=4, c(0.4, -0.6, 0.8, 0, 0.3, -0.7, -0.3, -1.2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#' 
#' 
#'  
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'  cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'  
#' # Set parameters of G0
#' 
#' priorMix <- list()
#' priorMix[["weights"]] <- p
#' priorMix[["parameters"]] <- list()  
#' for(j in 1:ncl){
#'  priorMix[["parameters"]][[j]] <- list()
#'  priorMix[["parameters"]][[j]][["b_xi"]] <- xi[,j]
#'  priorMix[["parameters"]][[j]][["b_psi"]] <- psi[,j]
#'  priorMix[["parameters"]][[j]][["lambda"]] <- crossprod(sdev[,,j])
#'  priorMix[["parameters"]][[j]][["kappa"]] <- 0.001
#'  priorMix[["parameters"]][[j]][["D_xi"]] <- 1
#'  priorMix[["parameters"]][[j]][["D_psi"]] <- 1
#'  priorMix[["parameters"]][[j]][["nu"]] <- d + 0.1
#' }
#'  
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'  
#'  # do some plots
#'  doPlot <- TRUE 
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
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'  MCMCsample_sn_EB <- DPMGibbsSkewN_EB(z, hyperG0=priorMix, a, b, N=1000, doPlot, plotevery=100, gg.add=list(theme_bw()))
#'  
#'  s <- summary(MCMCsample_sn, burnin = 500)
#'  print(s)
#'  plot(s)
#'  plot_ConvDPM(MCMCsample_sn, from=2)
#'  cluster_est_binder(MCMCsample_sn$c_list[50:500])
#'  
#'  
#'  
#'  
#'
#'
#'
#'
DPMGibbsSkewN_EB <- function (z, hyperG0, a, b, N, doPlot=TRUE, plotevery=1, diagVar=TRUE, ...){
    
    if(doPlot){library(ggplot2)}
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_xi <- matrix(0, nrow=p, ncol=n)
    U_psi <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    U_B = array(0, dim=c(2, 2, n))
    
    # U_SS is a list where each U_SS[[k]] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
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
    nbmix_prior <- length(hyperG0[["weights"]])
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    
    c <- sample(x=1:nbmix_prior, size=n, replace=TRUE)
    for (k in unique(c)){
        obs_k <- which(c==k)
        hyper_num <- k
        priormix <- hyperG0[["parameters"]][[hyper_num]]
        #cat("cluster ", k, ":\n")
        U_SS[[k]] <- update_SSsn(z=z[, obs_k], S=priormix, ltn=ltn[obs_k])
        NNiW <- rNNiW(U_SS[[k]], diagVar)
        U_xi[, k] <- NNiW[["xi"]]
        U_SS[[k]][["xi"]] <- NNiW[["xi"]]
        U_psi[, k] <- NNiW[["psi"]]
        U_SS[[k]][["psi"]] <- NNiW[["psi"]]
        U_Sigma[, , k] <- NNiW[["S"]]
        U_SS[[k]][["S"]] <- NNiW[["S"]]
        U_B[, ,k] <- U_SS[[k]][["B"]]
        m[k] <- length(obs_k)
    }
    
    
    
    if(is.null(hyperG0[["alpha"]])){
        alpha <- nbmix_prior/log(n)
    }else{
        alpha <- hyperG0[["alpha"]]
    }
    
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[1]] <- numeric(length(m))
    weights_list[[1]][unique(c)] <- table(c)/length(c)
    
    logposterior_list[[i]] <- 0#logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_B,
    #hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
    
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    
    if(doPlot){
        plot_DPMsn(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
    }else{
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    
    
    for(i in 2:N){
        nbClust <- length(unique(c))
        
        alpha <- c(alpha,
                   sample_alpha(alpha_old=alpha[i-1], n=n, 
                                K=nbClust, a=a, b=b)
        )
        slice <- sliceSampler_SkewN_EB(c=c, m=m, alpha=alpha[i], 
                                       z=z, hyperG0=hyperG0, 
                                       U_xi=U_xi, U_psi=U_psi, U_Sigma=U_Sigma,
                                       diagVar)
        m <- slice[["m"]]
        c <- slice[["c"]]        
        weights_list[[i]] <- slice[["weights"]]
        ltn <- slice[["latentTrunc"]]
        
        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in fullCl){
            obs_j <- which(c==j)
            #cat("cluster ", j, ":\n")
            
            #EB
            hyper_num <- sample(x=1:nbmix_prior, size=1, prob=hyperG0[[1]])
            priormix <- hyperG0[["parameters"]][[hyper_num]]
            
            U_SS[[j]] <- update_SSsn(z=z[, obs_j], S=priormix,  ltn=ltn[obs_j])
            NNiW <- rNNiW(U_SS[[j]], diagVar)
            U_xi[, j] <- NNiW[["xi"]]
            U_SS[[j]][["xi"]] <- NNiW[["xi"]]
            U_psi[, j] <- NNiW[["psi"]]
            U_SS[[j]][["psi"]] <- NNiW[["psi"]]
            U_Sigma[, , j] <- NNiW[["S"]]
            U_SS[[j]][["S"]] <- NNiW[["S"]]
            U_B[, ,j] <- U_SS[[j]][["B"]]
        }
        
        
        U_SS_list[[i]] <- U_SS[which(m!=0)]
        c_list[[i]] <- c
        
        logposterior_list[[i]] <- 0#logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_B,
        #hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
        
        cat(i, "/", N, " samplings:\n", sep="")
        cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
        
        if(doPlot && i/plotevery==floor(i/plotevery)){
            plot_DPMsn(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
        }else{
            cl2print <- unique(c)
            cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
        }
        
    }
    
    dpmclus <- list("mcmc_partitions" = c_list, 
                    "alpha"=alpha, 
                    "U_SS_list"=U_SS_list,
                    "weights_list"=weights_list, 
                    "logposterior_list"=logposterior_list, 
                    "data"=z,
                    "nb_mcmcit"=N,
                    "clust_distrib"="skewNormal")
    class(dpmclus) <- "DPMMclust"
    return(dpmclus)
}






