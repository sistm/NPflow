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
#'@param diagVar logical flag indicating wether the variance of each cluster is 
#'estimated as a diagonal matrix, or as a full matrix. 
#'Default is \code{TRUE} (diagonal variance).
#'
#'@param verbose logical flag indicating wether partition info is 
#'written in the console at each MCMC iteration.
#'
#'@author Boris Hejblum
#'
#'@export
#'
#'@examples
#'
#'
#'
#'
#' # Scaling up: ----
#' rm(list=ls())
#' #Number of data
#' n <- 2000
#' set.seed(1234)
#' 
#' # Sample data
#' d <- 10
#' nclust <- 20
#' m <- matrix(nrow=d, ncol=nclust, runif(d*nclust)*8)
#' # p: cluster probabilities
#' p <- runif(nclust)
#' p <- p/sum(p)
#'
#' # Covariance matrix of the clusters 
#' sdev <- array(dim=c(d, d, nclust))
#' for (j in 1:nclust){
#'     sdev[, ,j] <- matrix(NA, nrow=d, ncol=d)
#'     diag(sdev[, ,j]) <- abs(rnorm(n=d, mean=0.3, sd=0.1))
#'     sdev[, ,j][lower.tri(sdev[, ,j], diag = FALSE)] <- rnorm(n=d*(d-1)/2,
#'     mean=0, sd=0.05)
#'     sdev[, ,j][upper.tri(sdev[, ,j], diag = FALSE)] <- sdev[, ,j][lower.tri(sdev[, ,j], diag = FALSE)]
#' }
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'     c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'     z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'     cat(k, "/", n, " observations simulated\n", sep="")
#' }
#' 
#' # hyperprior on the Scale parameter of DPM
#' a <- 0.001
#' b <- 0.001
#' 
#' # Number of iterations
#' N <- 25
#' 
#' # do some plots
#' doPlot <- TRUE
#' 
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["mu"]] <- rep(0, d)
#' hyperG0[["kappa"]] <- 0.01
#' hyperG0[["nu"]] <- d + 2
#' hyperG0[["lambda"]] <- diag(d)/10
#' 
#' 
#' nbclust_init <- 30
#' MCMCsample <- gibbsDPMsliceprior(z, hyperG0, a, b, N, doPlot=F, nbclust_init=30)
#' 
#' plot_DPM(z=z, U_mu=MCMCsample$U_mu, U_Sigma=MCMCsample$U_Sigma,
#'          m=MCMCsample$partition, c=MCMCsample$clusters, 
#'          i=N, alpha=MCMCsample$alpha[[1000]], U_SS=MCMCsample$U_SS_list[[1000]],
#'          ellipses=T,
#'          dims2plot=c(1,2))
#' 
#'
DPMGibbsN_parallel <- function (Ncpus, type_connec,
                                z, hyperG0, a, b, N, doPlot=TRUE, 
                                nbclust_init=30, plotevery=1,
                                diagVar=TRUE, verbose=TRUE,
                                ...){
  
  if(doPlot){library(ggplot2)}
  library(doSNOW)
  
  p <- nrow(z)
  n <- ncol(z)
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
  
  #store log posterior probability
  logposterior_list <- list()
  
  
  m <- numeric(n) # number of obs in each clusters
  c <-numeric(n) # initial number of clusters
  
  # declare the cores
  cl <- makeCluster(Ncpus, type = type_connec)
  registerDoSNOW(cl)
  
  
  # Initialisation----
  # each observation is assigned to a different cluster
  # or to 1 of the 50 initial clusters if there are more than
  # 50 observations
  
  i <- 1
  if(n < nbclust_init){       
    c <- sample(x=1:ncol(z), size=n, replace=FALSE)
  }else{
    c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
  }
  
  browser()
  
  res <- foreach(k=unique(c), .packages="NPflow")%dopar%{
    obs_k <- which(c==k)
    U_SS_par <- update_SS(z=z[, obs_k], S=hyperG0)
    NiW <- rNiW(U_SS_par, diagVar)
    U_mu_par <- NiW[["mu"]]
    U_Sigma_par <- NiW[["S"]]
    m_par <- length(obs_k)
    list("U_SS_par"=U_SS_par, "U_mu_par"=U_mu_par, 
         "U_Sigma_par"=U_Sigma_par, "m_par"=m_par)
  }
  
  U_SS <- lapply(res, FUN="[[", "U_SS_par")
  names(U_SS) <- unique(c)
  U_mu <- lapply(res, FUN="[[", "U_mu_par")
  names(U_mu) <- unique(c)
  U_Sigma <- lapply(res, FUN="[[", "U_Sigma_par")
  names(U_Sigma) <- unique(c)
  m[unique(c)] <- unlist(lapply(res, FUN="[[", "m_par"))
  
  
  
  alpha <- log(n)
  
  U_SS_list[[i]] <- U_SS
  c_list[[i]] <- c
  weights_list[[1]] <- numeric(length(m))
  weights_list[[1]][unique(c)] <- table(c)/length(c)
  
  logposterior_list[[i]] <- logposterior_DPMG(z=z, mu=U_mu, Sigma=U_Sigma,
                                              hyper=hyperG0, m=m, 
                                              c=c, alpha=alpha[i], 
                                              n=n, a=a, b=b)
  
  if(verbose){    
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    cl2print <- unique(c)
    cat("  ",length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
  }
  
  if(doPlot){
    plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, 
             m=m, c=c, i=i, alpha=alpha[length(alpha)], U_SS=U_SS, ...)
  }
  
  
  
  
  # MCMC iterations
  
  for(i in 2:N){
    nbClust <- length(unique(c))
    alpha <- c(alpha,
               sample_alpha(alpha_old=alpha[i-1], n=n, 
                            K=nbClust, a=a, b=b)
    )
    slice <- slice_sample_parallel(c=c, m=m, alpha=alpha[i], 
                                   z=z, hyperG0=hyperG0, 
                                   U_mu=U_mu, U_Sigma=U_Sigma, diagVar=diagVar)
    m <- slice[["m"]]
    c <- slice[["c"]]
    weights_list[[i]] <- slice[["weights"]]
    
    
    # Update cluster locations
    fullCl <- which(m!=0)
    res <- foreach(k=fullCl, .packages="NPflow")%dopar%{
      obs_k <- which(c==k)
      U_SS_par <- update_SS(z=z[, obs_k], S=hyperG0)
      NiW <- rNiW(U_SS_par, diagVar)
      U_mu_par <- NiW[["mu"]]
      U_Sigma_par <- NiW[["S"]]
      list("U_SS_par"=U_SS_par, "U_mu_par"=U_mu_par, 
           "U_Sigma_par"=U_Sigma_par)
    }
    
    U_SS <- lapply(res, FUN="[[", "U_SS_par")
    names(U_SS) <- as.character(fullCl)
    U_mu <- lapply(res, FUN="[[", "U_mu_par")
    names(U_mu) <- as.character(fullCl)
    U_Sigma <- lapply(res, FUN="[[", "U_Sigma_par")
    names(U_Sigma) <- as.character(fullCl)
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    
    
    logposterior_list[[i]]  <- logposterior_DPMG(z=z, mu=U_mu, Sigma=U_Sigma,
                                                 hyper=hyperG0, m=m, 
                                                 c=c, alpha=alpha[i], 
                                                 n=n, a=a, b=b)
    
    if(verbose){    
      cat(i, "/", N, " samplings:\n", sep="")
      cat("  logposterior = ", sum(logposterior_list[[i]]) , "\n", sep="")
      cl2print <- unique(c)
      cat("  ",length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    if(doPlot && i/plotevery==floor(i/plotevery)){
      plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, 
               m=m, c=c, i=i, alpha=alpha[length(alpha)], U_SS=U_SS, ...)
    }
    
  }
  
  stopCluster(cl)
  
  return(list("clusters" = c, "U_mu" = U_mu, "U_Sigma" = U_Sigma, 
              "partition" = m, "alpha"=alpha, "U_SS_list"=U_SS_list,
              "c_list" = c_list, "weights_list"=weights_list, 
              "logposterior_list"=logposterior_list,
              "nb_mcmcit"=N,
              "clust_distrib"="Normal",
              "hyperG0"=hyperG0))
}






