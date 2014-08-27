#'Slice Sampling of Dirichlet Process Mixture of skew  Students's t-distibutions
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
#'@param verbose logical flag indicating wether partition info is 
#'written in the console at each MCMC iteration.
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
#' #set.seed(4321)
#' 
#' 
#' d <- 2
#' ncl <- 4
#' 
#' # Sample data
#' 
#' 
#' xi <- matrix(nrow=d, ncol=1, c(-2.5, -3))
#' psi <- matrix(nrow=d, ncol=1, c(-0.3, -0.7))
#' nu <- 27
#' sdev <- .5*diag(2)
#' 
#'  
#' c <- rep(1,n)
#' w <- rep(1,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  w <- rgamma(n, shape=nu/2, rate=nu/2)
#'  z[,k] <- xi + psi*rtruncnorm(n=1, a=0, b=Inf, mean=0, sd=1/sqrt(w[k])) + (sdev/sqrt(w[k]))%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'  cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'  
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rowMeans(z)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d+1
#' hyperG0[["lambda"]] <- diag(apply(z,MARGIN=1, FUN=var))
#'  
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'  
#'  # do some plots
#'  doPlot <- FALSE 
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
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'  MCMCsample_st <- DPMGibbsSkewT_test1K(z, hyperG0, a, b, N=1000, 
#'  doPlot, nbclust_init, diagVar=FALSE)
#'  sapply(MCMCsample_st$U_SS_list, "[[", "df")
#'  sapply(MCMCsample_st$U_SS_list, "[[", "xi")
#'  sapply(MCMCsample_st$U_SS_list, "[[", "psi")
#'  plot(density(sapply(MCMCsample_st$U_SS_list, "[[", "df")[500:1000]))
#'  plot(density(sapply(MCMCsample_st$U_SS_list, "[[", "psi")[2,500:1000]))


DPMGibbsSkewT_test1K <- function (z, hyperG0, a, b, N, doPlot=TRUE, 
                                  nbclust_init=30, plotevery=1, 
                                  diagVar=TRUE, verbose=TRUE,
                                  ...){
    
    if(doPlot){library(ggplot2)}
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_xi <- matrix(0, nrow=p, ncol=n)
    U_psi <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    U_df = 10
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
    
    
    m <- numeric(n) # number of obs in each clusters
    c <- numeric(n) # cluster label of ech observation
    ltn <- rtruncnorm(n, a=0, b=Inf, mean=0, sd=1) # latent truncated normal
    sc <- rep(1,n)
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    
    U_SS <- update_SSst(z=z, S=hyperG0, ltn=ltn, scale=sc, df=U_df)
    NNiW <- rNNiW(U_SS, diagVar)
    U_xi <- NNiW[["xi"]]
    U_SS[["xi"]] <- NNiW[["xi"]]
    U_psi <- NNiW[["psi"]]
    U_SS[["psi"]] <- NNiW[["psi"]]
    U_Sigma <- NNiW[["S"]]
    U_SS[["S"]] <- NNiW[["S"]]
    U_B<- U_SS[["B"]]
    
    alpha <- c(log(n))
    
    
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[1]] <- numeric(length(m))
    weights_list[[1]][unique(c)] <- table(c)/length(c)
    
    if(verbose){
        cat(i, "/", N, " samplings:\n", sep="")
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    acc_rate <- 0
    
    if(N>1){
        for(i in 2:N){
            nbClust <- length(unique(c))
            
            alpha <- c(alpha,
                       sample_alpha(alpha_old=alpha[i-1], n=n, 
                                    K=nbClust, a=a, b=b)
            )
            
            slice <- sliceSampler_SkewT_test1K(c=c, m=m, alpha=alpha[i], 
                                               z=z, hyperG0=hyperG0, 
                                               U_xi=U_xi, U_psi=U_psi, 
                                               U_Sigma=U_Sigma, U_df=U_df,
                                               scale=sc, diagVar)
            
            ltn <- slice[["latentTrunc"]]
            
            
            # Update cluster locations            
            U_SS <- update_SSst(z=z, S=hyperG0, 
                                ltn=ltn, scale=sc, 
                                df=U_df)
            NNiW <- rNNiW(U_SS, diagVar)
            U_xi<- NNiW[["xi"]]
            U_SS[["xi"]] <- NNiW[["xi"]]
            U_psi <- NNiW[["psi"]]
            U_SS[["psi"]] <- NNiW[["psi"]]
            U_Sigma <- NNiW[["S"]]
            U_SS[["S"]] <- NNiW[["S"]]
            U_B <- U_SS[["B"]]
            
            
            update_scale <- sample_scale_test1K(c=c, m=m, z=z, U_xi=U_xi, 
                                                U_psi=U_psi, U_Sigma=U_Sigma, 
                                                U_df=U_df, ltn=ltn, 
                                                weights=weights_list[[i]],
                                                scale=sc)
            U_df <- update_scale[["df"]]
            sc <- update_scale[["scale"]]
            acc_rate <- acc_rate + update_scale[["acc_rate"]]
            U_SS[["df"]] <- U_df
            
            
            U_SS_list[[i]] <- U_SS
            c_list[[i]] <- c
            
            if(verbose){
                cat(i, "/", N, " samplings:\n", sep="")
                cl2print <- unique(c)
                cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
            }
        }
    }
    acc_rate <- acc_rate/N
    
    dpmclus <- list("mcmc_partitions" = c_list, 
                    "alpha"=alpha, 
                    "U_SS_list"=U_SS_list,
                    "weights_list"=weights_list, 
                    "data"=z,
                    "nb_mcmcit"=N,
                    "clust_distrib"="skewT",
                    "acc_rate"=acc_rate)
    class(dpmclus) <- "DPMMclust"
    return(dpmclus)
}






