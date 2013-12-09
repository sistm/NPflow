#'Gibbs Sampling with Alghorithm 2
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
#'@export gibbsDPMalgo2
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
#'  GibSample <- gibbsDPMalgo2(z, hyperG0, alpha, N, doPlot)
#'
#'
gibbsDPMalgo2 <- function (z, hyperG0, alpha, N, doPlot=TRUE){
    
    if(doPlot){library(ggplot2)}
    
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    
    # U_SS is a list where each U_SS[k] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    m <- numeric(n)
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
        plot_DPM2(z, U_mu, m, c, i)
    }
    
    
    for(i in 2:N){
        for (k in 1:n){
            # Update cluster assignments c
            
            m[c[k]] <- m[c[k]] - 1
            U_SS[[c[k]]] <- downdate_SS(z[, k], U_SS[[c[k]]])
            c[k] <- sample_c(m, alpha, z[, k], hyperG0, U_mu, U_Sigma)
            m[c[k]] <- m[c[k]] + 1
            
            if(m[c[k]]>1){
                U_SS[[c[k]]] <- update_SS(z[, k], U_SS[[c[k]]])
            } else {
                U_SS[[c[k]]] <- update_SS(z[, k], hyperG0)
                NiW <- normalinvwishrnd(U_SS[[c[k]]])
                U_mu[, c[k]] <- NiW[["mu"]]
                U_Sigma[, , c[k]] <- NiW[["S"]]
            }
        }
        
        # Update cluster locations U
        ind <- which(m!=0)
        
        for(j in 1:length(ind)){
            NiW <- normalinvwishrnd(U_SS[[ind[j]]])
            U_mu[, ind[j]] <- NiW[["mu"]]
            U_Sigma[, , ind[j]] <- NiW[["S"]]
        }
        
        cat(i, "/", N, " samplings\n", sep="")
        if(doPlot){
            plot_DPM2(z, U_mu, m, c, i)
        }
    }
    return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
                "partition" = m))
}








# Subfunctions ----
sample_c <- function(m, alpha, z, hyperG0, U_mu, U_Sigma){
    
    fullCl <- which(m!=0) # indexes of non empty clusters
    r <- sum(m)
    n <- numeric(length(fullCl))
    for (i in 1:length(fullCl)){
        n[i] <- mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , 
                       mean = U_mu[, fullCl[i]], 
                       varcovM = U_Sigma[, , fullCl[i]])*m[fullCl[i]]  
    }
    
    n0 <- pred(z, hyperG0)
    const <- sum(n)/(alpha + r) + alpha/(alpha + r)*n0
    p0 <- alpha/(alpha + r)*n0/const # probability of sampling a new item
    
    u <- runif(n=1, min = 0, max = 1)
    if (u<p0){
        # Accept: allocate to a new cluster
        # cat("acceptation:", u, "<", p0, "\n")
        K <- which(m==0)[1]
    } 
    else{
        # Reject: allocate to a previously non empty cluster
        # cat("rejection:", u, ">=", p0, "\n")
        u1  <-  u - p0
        ind <- which(cumsum(n/const/(alpha+r))>rep(u1,length(n)))[1]
        K <- fullCl[ind]
    }
    return(K)
}


plot_DPM2 <- function(z, U_mu, m, c, i){
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
