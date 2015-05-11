#'Slice Sampling of the Dirichlet Process Mixture Model
#'with a prior on alpha
#'
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
#'@export 
#'
#'@examples
#' rm(list=ls())
#' library(ggplot2)
#' #Number of data
#' n <- 500
#' #n <- 2000
#' set.seed(1234)
#' #set.seed(123)
#' #set.seed(4321)
#' 
#' # Sample data
#' m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' 
#' sdev <- array(dim=c(2,2,4))
#' sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#' c <- rep(0,n)
#' z <- matrix(0, nrow=2, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'  
#'  d<-2
#'  # Set parameters of G0
#'  hyperG0 <- list()
#'  hyperG0[["mu"]] <- rep(0,d)
#'  hyperG0[["kappa"]] <- 0.001
#'  hyperG0[["nu"]] <- d+2
#'  hyperG0[["lambda"]] <- diag(d)/10
#'  
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'  
#'  # Number of iterations
#'  N <- 30
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
#'        + ggtitle("Toy example Data"))
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
#'  Rprof("Rprof.out")
#'  MCMCsample <- DPMGibbsN(z, hyperG0, a, b, N=500, doPlot, nbclust_init, plotevery=50)
#'                  
#'  plot_ConvDPM(MCMCsample, from=2)
#'  
#'  s <- summary(MCMCsample, burnin = 10, thin=2, posterior_approx=FALSE, K=0,dist="Normal",
#'  lossFn = "MBinderN",lambda=0,a=1,b=1)
#'  
#'  F <- FmeasureC(pred=s$point_estim$c_est, ref=c)
#'  
#'  
#'  # k-means
#'  
#'  plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster,
#'       xlab = "d = 1", ylab= "d = 2", main="k-means with K=4 clusters")
#'       
#'  KM <- kmeans(t(z), centers=4)
#'  dataKM <- data.frame("X"=z[1,], "Y"=z[2,], 
#'                     "Cluster"=as.character(KM$cluster))
#'  dataCenters <- data.frame("X"=KM$centers[,1], 
#'                            "Y"=KM$centers[,2], 
#'                            "Cluster"=rownames(KM$centers))
#'  
#'  p <- (ggplot(dataKM) 
#'        + geom_point(aes(x=X, y=Y, col=Cluster))
#'        + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster), 
#'                     data=dataCenters, shape=22, size=5)
#'        + scale_colour_discrete(name="Cluster")
#'        + ggtitle("K-means with K=4 clusters\n"))
#'  p
#'  
#'  postalpha <- data.frame("alpha"=MCMCsample$alpha[50:500], 
#'                          "distribution" = factor(rep("posterior",500-49), 
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(postalpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..), binwidth=.1,
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of alpha\n")
#'        + geom_vline(aes(xintercept=mean(alpha, na.rm=T)),   # Ignore NA values for mean
#'                     color="red", linetype="dashed", size=1)  # Overlay with transparent density plot            
#'      )
#'  p
#'  
#'  p <- (ggplot(drop=FALSE, alpha=.6)
#'        + geom_density(aes(x=alpha, fill=distribution), 
#'                       color=NA, alpha=.6,
#'                       data=prioralpha)
#'        + geom_density(aes(x=alpha, fill=distribution), 
#'                       color=NA, alpha=.6,
#'                       data=postalpha)
#'        + ggtitle("Prior and posterior distributions of alpha\n")
#'        + scale_fill_discrete(drop=FALSE)
#'      )
#'  p
#'
#'
#'
#'
#'

DPMGibbsN <- function (z, hyperG0, a, b, N, doPlot=TRUE, 
                                nbclust_init=30, plotevery=1, ...){
    
    if(doPlot){library(ggplot2)}
    
    p <- nrow(z)
    n <- ncol(z)
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    listU_mu<-list()
    listU_Sigma<-list()
    
    # U_SS is a list where each U_SS[[k]] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    #store U_SS :
    U_SS_list <- list()
    #store clustering :
    c_list <- list()
    #store sliced weights
    weights_list <- list()
    
    #store concentration parameter
    #alpha <- list()
    
    #store log posterior probability
    logposterior_list <- list()
    
    m <- numeric(n) # number of obs in each clusters
    c <-numeric(n)
    # initial number of clusters
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    if(n<nbclust_init){       
        for (k in 1:n){
            c[k] <- k
            #cat("cluster ", k, ":\n")
            U_SS[[k]] <- update_SS(z=z[, k, drop=FALSE], S=hyperG0)
            NiW <- rNiW(U_SS[[k]],diagVar=FALSE)
            
            U_mu[, k] <- NiW[["mu"]]
            U_SS[[k]][["mu"]] <- NiW[["mu"]]
            
            U_Sigma[, , k] <- NiW[["S"]]
            U_SS[[k]][["S"]] <- NiW[["S"]]
            
            m[k] <- m[k]+1
            U_SS[[k]][["weight"]] <- 1/n
        }
    } else{
        c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
        for (k in unique(c)){
            obs_k <- which(c==k)
            #cat("cluster ", k, ":\n")
            U_SS[[k]] <- update_SS(z=z[, obs_k,drop=FALSE], S=hyperG0)
            NiW <- rNiW(U_SS[[k]],diagVar=FALSE)
            
            U_mu[, k] <- NiW[["mu"]]
            U_SS[[k]][["mu"]] <- NiW[["mu"]]
            
            U_Sigma[, , k] <- NiW[["S"]]
            U_SS[[k]][["S"]] <- NiW[["S"]]
            
            m[k] <- length(obs_k)
            U_SS[[k]][["weight"]] <- m[k]/n
        }
    }
    listU_mu[[i]]<-U_mu
    listU_Sigma[[i]]<-U_Sigma
    
    alpha <- c(log(n))
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[i]] <- numeric(length(m))
    weights_list[[i]][unique(c)] <- table(c)/length(c)
    
    logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma, 
                                                hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
    
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    
    if(doPlot){
        plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, 
                 m=m, c=c, i=i, alpha=alpha[[i]], U_SS=U_SS, ...)
    }else{
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
   
    for(i in 2:N){
        nbClust <- length(unique(c))
        alpha <- c(alpha,sample_alpha(alpha_old=alpha[i-1], n=n, 
                                K=nbClust, a=a, b=b))
        
        slice <- slice_sample(c=c, m=m, alpha=alpha[i], 
                              z=z, hyperG0=hyperG0, 
                              U_mu=U_mu, U_Sigma=U_Sigma)
        m <- slice[["m"]]
        c <- slice[["c"]]        
        weights_list[[i]] <- slice[["weights"]]
        U_mu<-slice[["U_mu"]]
        U_Sigma<-slice[["U_Sigma"]]

        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in fullCl){
            obs_j <- which(c==j)
            #cat("cluster ", j, ":\n")
            U_SS[[j]] <- update_SS(z=z[, obs_j,drop=FALSE], S=hyperG0)
            NiW <- rNiW(U_SS[[j]],diagVar=FALSE)
            
            U_mu[, j] <- NiW[["mu"]]
            U_SS[[j]][["mu"]] <- NiW[["mu"]]
            
            U_Sigma[, , j] <- NiW[["S"]]
            U_SS[[j]][["S"]] <- NiW[["S"]]
            
            U_SS[[j]][["weight"]] <- weights_list[[i]][j]
            #cat("sampled S =", NiW[["S"]], "\n\n\n")
        }
        
        listU_mu[[i]]<-U_mu
        listU_Sigma[[i]]<-U_Sigma
        U_SS_list[[i]] <- U_SS[which(m!=0)]
        c_list[[i]] <- c
        logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma, 
                                                    hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)
        
        cat(i, "/", N, " samplings:\n", sep="")
        cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
        
        if(doPlot && i/plotevery==floor(i/plotevery)){
            plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, m=m, c=c, i=i,
                     alpha=alpha[i], U_SS=U_SS, ...)
        }else{
            cl2print <- unique(c)
            cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
        }
        
    }
    
#     return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
#                 "partition" = m, "alpha"=alpha, "U_SS_list"=U_SS_list,
#                 "c_list" = c_list, "weights_list"=weights_list, 
#                 "logposterior_list"=logposterior_list,
#                 "nb_mcmcit"=N,
#                 "clust_distrib"="Normal",
#                 "hyperG0"=hyperG0))
dpmclus <- list("mcmc_partitions" = c_list, 
                "alpha"=alpha, 
#                 "U_mu" = U_mu, 
#                 "U_Sigma" = U_Sigma,
                "listU_mu"=listU_mu,
                "listU_Sigma"=listU_Sigma,
                "U_SS_list"=U_SS_list,
                "weights_list"=weights_list, 
                "logposterior_list"=logposterior_list, 
                "data"=z,
                "nb_mcmcit"=N,
                "clust_distrib"="Normal",
                "hyperG0"=hyperG0)
class(dpmclus) <- "DPMMclust"
return(dpmclus)
}






