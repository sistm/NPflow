#'Burning mcmc iterations from a Dirichlet Process Mixture Model
#'
#'@name burn.DPMMclust
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin the number of MCMC iterations to burn (defaults is half).
#'
#'@param gs optionnal vector of length \code{n} containing the gold standard 
#'partition of the \code{n} observations to compare to the point estimate.
#'
#'@param ... further arguments passed to or from other methods.
#'
#'@return a \code{DPMMclust} object minus the burnt iterations
#'
#'@author Boris Hejblum
#'
#'@export burn.DPMMclust
#'
#'@seealso \link{summary.DPMMclust}
#'

burn.DPMMclust <- function(x, burnin=0){
    
    if(burnin>=length(x[["mcmc_partitions"]])){
        stop("burnin argument is superior to the number of MCMC iterations sampled")
    }
    
    xburnt <- list()
    
    N <- x[["nb_mcmcit"]]
    
    xburnt[["mcmc_partitions"]] <- x[["mcmc_partitions"]][(burnin+1):N]
    xburnt[["alpha"]] <- x[["alpha"]][(burnin+1):N]
    xburnt[["U_SS_list"]] <- x[["U_SS_list"]][(burnin+1):N]
    xburnt[["weights_list"]] <- x[["weights_list"]][(burnin+1):N]
    xburnt[["logposterior_list"]] <- x[["logposterior_list"]][(burnin+1):N]
    xburnt[["data"]] <- x[["data"]]
    xburnt[["nb_mcmcit"]] <- N-burnin
    xburnt[["clust_distrib"]] <- x[["clust_distrib"]]
    
    class(xburnt) <- "DPMMclust"
    
    return(xburnt)
    
}