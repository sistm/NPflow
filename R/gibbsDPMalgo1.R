#'Gibbs Sampling with Alghorithm 1
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
#'@export gibbsDPMalgo1
#'
#'@examples
#' rm(list=ls())
#' #Number of data
#' n <- 100
#' set.seed(1231)
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
#'  GibSample <- gibbsDPMalgo1(z, hyperG0, alpha, N, doPlot)
#'
#'
gibbsDPMalgo1 <- function (z, hyperG0, alpha, N, doPlot=TRUE){
    
    if(doPlot){library(ggplot2)}

    
    p <- dim(z)[1]
    n <- dim(z)[2]
    theta_mu <- matrix(0, nrow=p, ncol=n)
    theta_Sigma = array(0, dim=c(p, p, n))
    
    
    # Initialisation----
    i=1
    hyper <- update_SS(z[, 1], hyperG0)
    NiW <- normalinvwishrnd(hyper)
    theta_mu[, 1] <- NiW[["mu"]]
    theta_Sigma[, , 1] <- NiW[["S"]]
    for (k in 2:n){
        sampTheta <- sample_theta(alpha, z=z[, k], hyperG0=hyperG0, 
                                  theta_mu_notk=theta_mu[, 1:(k-1)], 
                                  theta_Sigma_notk=theta_Sigma[, , 1:(k-1)])
        theta_mu[, k] <- sampTheta[["mu"]]
        theta_Sigma[, , k] = sampTheta[["Sigma"]]
    }
    
    cat(i, "/", N, " samplings\n", sep="")
    if(doPlot){
        plot_DPM1(z, theta_mu, n, i)
    }
    
    
    for(i in 2:N){
        for (k in 1:n){
            # Sample theta_k | theta_{-k}, z_k
            #ind_notk <- c(1:n)[-k]
            sampTheta <- sample_theta(alpha, z=z[, k], hyperG0, 
                                      theta_mu_notk=theta_mu[, -k], 
                                      theta_Sigma_notk=theta_Sigma[, , -k])
            theta_mu[, k] <- sampTheta[["mu"]]
            theta_Sigma[, , k] = sampTheta[["Sigma"]]
        }
        
        cat(i, "/", N, " samplings\n", sep="")
        if(doPlot){
            plot_DPM1(z, theta_mu, n, i)
        }
    }
    return(list("mu"=theta_mu, "Sigma"=theta_Sigma))
}








# Subfunctions ----

sample_theta <- function(alpha, z, hyperG0, theta_mu_notk, theta_Sigma_notk){

    if(is.null(dim(theta_mu_notk))){
        p <- 1
        n <-mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , mean = theta_mu_notk, varcovM = theta_Sigma_notk)
    }else{
        p <- dim(theta_mu_notk)[2]
        n <- numeric(p)
        for (i in 1:p){
            n[i] <- mvnpdf(x = matrix(z, nrow= 1, ncol=length(z)) , mean = theta_mu_notk[,i], varcovM = theta_Sigma_notk[, , i])  
        }
    }
    
    n0 <- pred(z, hyperG0)
    const <- alpha*n0+sum(n)
    p0 <- alpha*n0/const
    #cat("probability of allocating to a new cluster:", p0, "\n")

    u <- runif(n=1, min = 0, max = 1)
    if (u<p0){
        # Accept: sample new value
        # cat("acceptation:", u, "<", p0, "\n")
        hyper <- update_SS(z, hyperG0)
        NiW <- normalinvwishrnd(hyper)
        theta_mu_k <- NiW[["mu"]]
        theta_Sigma_k <- NiW[["S"]]
    } 
    else{
        # Reject: sample old value
        # cat("rejection:", u, ">=", p0, "\n")
        u1  <-  u - p0
        ind <- which(cumsum(n/const)>rep(u1,length(n)))[1]
        
        if(is.null(dim(theta_mu_notk))){
            theta_mu_k = theta_mu_notk  
            theta_Sigma_k = theta_Sigma_notk          
        } else{
            theta_mu_k = theta_mu_notk[,ind]
            theta_Sigma_k = theta_Sigma_notk[, , ind] 
        }  
    }
    return(list("mu"=theta_mu_k, "Sigma"= theta_Sigma_k))
}

#calculer le log de la posterior à chaque iteration
#doit augmenter puis stagner
#MAP =>estimaton des param pour le max


plot_DPM1 <- function(z, theta_mu, n, i){
    ind <- unique(theta_mu[1,])
    cl <- numeric(n)
    U_mu <- matrix(0, nrow=2, ncol=length(ind))
    for(j in 1:length(ind)){
        ind2 <- which(theta_mu[1, ]==ind[j])
        cl[ind2] <- j
        U_mu[, j] <- theta_mu[, ind2[1]]
    }
    z2plot <- cbind.data.frame("X"=z[1,],"Y"=z[2,],"Cluster"=as.factor(cl))
    U2plot <- cbind.data.frame("X"=U_mu[1,],"Y"=U_mu[2,],"Cluster"=factor(1:dim(U_mu)[2]))
    p <- (ggplot(z2plot) 
          + geom_point(aes(x=X, y=Y, col=Cluster), data=z2plot) 
          + geom_point(aes(x=X, y=Y, col=Cluster), data=U2plot, shape="X", size=5)
          + ggtitle(paste("Gibbs sampling for DPM - algo 1\nIteration", i))
    )
    print(p)
}
