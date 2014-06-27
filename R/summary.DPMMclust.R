#'Summarizing Dirichlet Process Mixture Models
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
#'@return a \code{list}: 
#'  \itemize{
#'      \item{\code{c_est}:}{ a vector of length \code{n}. Point estimate of the partition}
#'      \item{\code{cost}:}{ a vector of length \code{N}. \code{cost[j]} is the cost 
#' associated to partition \code{c[[j]]}}
#'      \item{\code{similarity}:}{  matrix of size \code{n x n}. Similarity matrix 
#' (see \link{similarityMat})}
#'  }
#'
#'@author Boris Hejblum
#'
#'@export summary.DPMMclust
#'
#'@export print.summaryDPMMclust
#'
#'@export plot.summaryDPMMclust
#'
#'@seealso \link{similarityMat}
#'

summary.DPMMclust <- function(x, burnin=0, gs=NULL,...){
    
    x_invar <- burn.DPMMclust(x, burnin = burnin)
    point_estim <- cluster_est_binder(x_invar$mcmc_partitions)
    index_estim <- point_estim$opt_ind
    loss <- NA
    if(!is.null(gs)){
        loss <- evalClust_Binder(c=point_estim$c_est, gs=gs, ...)
    }
    
    s <- c(x_invar, list("burnin"=burnin, 
                   "point_estim"=point_estim,
                   "loss"=loss,
                   "index_estim"=index_estim))
    class(s) <- "summaryDPMMclust"
    
    invisible(s)
}

print.summaryDPMMclust <- function(s,...){
    
    cat(class(s), "object with", s$nb_mcmcit, "MCMC iterations:\n")
    cat(rep("-",40),"\n", sep="")
    cat(rep("-",40),"\n", sep="")
    cat("Burnin =", s$burnin, "iterations\n\n")
    cat("Point estimate of the partition with a ", s$clust_distrib," mixture:\n", sep="")
    t<-table(s$point_estim$c_est)
    t2print <- paste(formatC(as.vector(t/length(s$point_estim$c_est)), 
                             digits=2),"%", sep="")
    names(t2print) <- names(t)
    print(t2print, quote=F)
    cat("\nLoss of the point estimate partition:", s$loss, "\n")
    cat(rep("-",40),"\n", sep="")
    
}

plot.summaryDPMMclust <- function(s, ...){
    
    plot_ConvDPM(s, shift=s$burnin)
        
    ind <- s$index_estim
    
    if(s$clust_distrib=="Normal"){
        
    }else if(s$clust_distrib=="skewNormal"){
        cat("This plotting may take a few sec... ")
        plot_DPMsn(z=s$data,
                   c=s$point_estim$c_est, 
                   i=ind, 
                   alpha=s$alpha[ind], 
                   U_SS=s$U_SS_list[[ind]], 
                   ellipses=TRUE,
                   gg.add=list(theme_bw()),
                   nbsim_dens=200000,
                   ...
        )
        cat("DONE!")
    }
    
}
