#'Slice Sampling of the Dirichlet Process Mixture Model
#'with a prior on alpha
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
#'@export gibbsDPMsliceprior
#'
#'@examples
#' rm(list=ls())
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
#'  a <- 1
#'  b <- 5
#'  prioralpha <- data.frame("alpha"=rgamma(n=5000, a, 1/b), 
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
#'  # Number of iterations
#'  N <- 30
#'  
#'  # do some plots
#'  doPlot <- TRUE 
#'  nbclust_init <- 30
#'  
#'  p <- (ggplot(data, aes(x=X, y=Y)) 
#'        + geom_point()
#'        + ggtitle("Toy example Data"))
#'  p
#'  
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  MCMCsample <- gibbsDPMsliceprior(z, hyperG0, a, b, N, doPlot, nbclust_init)
#'  
#'  plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster,
#'       xlab = "d = 1", ylab= "d = 2", main="k-means with K=4 clusters")
#'       
#'  KM <- kmeans(t(z), centers=4)
#'  data <- data.frame("X"=z[1,], "Y"=z[2,], 
#'                     "Cluster"=as.character(KM$cluster))
#'  dataCenters <- data.frame("X"=KM$centers[,1], 
#'                            "Y"=KM$centers[,2], 
#'                            "Cluster"=rownames(KM$centers))
#'  
#'  p <- (ggplot(data) 
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
#' # Scaling up: ----
#' #'rm(list=ls())
#' #Number of data
#' n <- 100000
#' #n <- 50000
#' set.seed(1234)
#' #set.seed(123)
#' #set.seed(4321)
#'
#' # Sample data
#' d <- 15
#' nclust <- 80
#' m <- matrix(nrow=d, ncol=nclust, runif(d*nclust)*8)
#' # p: cluster probabilities
#' p <- runif(nclust)
#' p <- p/sum(p)
#'
#'
#'
#' sd <- array(dim=c(d, d, nclust))
#' for (j in 1:nclust){
#'     sd1 <- runif(1)*0.5
#'     sd2 <- runif(1)*0.05
#'     sd[, ,j] <- matrix(nrow=d, ncol=d, c(sd1, sd2, sd2, sd1))
#' }
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'     c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'     z[,k] <- m[, c[k]] + sd[, , c[k]]%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'     cat(k, "/", n, " observations simulated\n", sep="")
#' }
#' 
#' plot(z[1,], z[2,])
#' 
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["mu"]] <- rep(0, d)
#' hyperG0[["kappa"]] <- 1
#' hyperG0[["nu"]] <- nclust
#' hyperG0[["lambda"]] <- diag(d)
#' 
#' # hyperprior on the Scale parameter of DPM
#' a <- 0.1
#' b <- 0.001
#' plot(density(rgamma(n=5000, a, 1/b)))
#' # Number of iterations
#' N <- 25
#' 
#' # do some plots
#' doPlot <- TRUE
#' nbclust_init <- 200
#' 
#' MCMCsample <- gibbsDPMsliceprior(z, hyperG0, a, b, N, doPlot, nbclust_init=100)
#' 
#'
gibbsDPMsliceprior <- function (z, hyperG0, a, b, N, doPlot=TRUE, nbclust_init=30){
    
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
    #store sliced weights
    weights_list <- list()
    
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
    
    
    
    alpha <- log(n)
    
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[1]] <- numeric(length(m))
    weights_list[[1]][1:length(unique(c))] <- table(c)/length(c)
    
    
    cat(i, "/", N, " samplings\n", sep="")
    if(doPlot){
        plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, 
                 m=m, c=c, i=i, alpha=alpha[length(alpha)])
    }
    
    
    for(i in 2:N){
        
        nbClust <- length(unique(c))
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
        
        weights_list[[i]] <- slice[["weights"]]
        
        
        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in fullCl){
            obs_j <- which(c==j)
            U_SS[[j]] <- update_SS(z=z[, obs_j], S=hyperG0)
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
                "c_list" = c_list, "weights_list"=weights_list))
}






