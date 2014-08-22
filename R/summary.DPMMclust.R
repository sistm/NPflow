#'Summarizing Dirichlet Process Mixture Models
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin the number of MCMC iterations to burn (defaults is half)
#'
#'@param gs optionnal vector of length \code{n} containing the gold standard 
#'partition of the \code{n} observations to compare to the point estimate
#'
#'@param ... further arguments passed to or from other methods
#'
#'@return a \code{list}: 
#'  \itemize{
#'      \item{\code{burnin}:}{an integer passing along the \code{burnin} argument}
#'      \item{\code{point_estim}:}{}
#'      \item{\code{loss}:}{}
#'      \item{\code{index_estim}:}{}
#'  }
#'
#'@author Boris Hejblum
#'
#'@export 
#'
#'@importFrom gplots heatmap.2
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


#' Methods for a summary of a 'DPMMclust' object
#'@rdname methods.summaryDPMMclust
#'@aliases print.summaryDPMMclust, plot.summaryDPMMclust
#'@param s a \code{summaryDPMMclust} object
#'@param hm logical flag to plot the heatmap of the similarity matrix. 
#'Default is \code{FALSE}.
#'@param nbsim_densities the number of simulated observations to be used
#'to plot the density lines of the clusters in the point estimate partition plot
#'@param ... further arguments passed to or from other methods
#'@author Boris Hejblum
#'@export
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

#'@export
#'@rdname methods.summaryDPMMclust
plot.summaryDPMMclust <- function(s, hm=FALSE, nbsim_densities=50000, ...){
    plot_ConvDPM(s, shift=s$burnin)
    
    ind <- s$index_estim
    
    
    cat("Plotting point estimate (may take a few sec)... ")
    if(s$clust_distrib=="Normal"){
        plot_DPM(z=s$data,
                 c=s$point_estim$c_est, 
                 i=ind+s$burnin, 
                 alpha=s$alpha[ind], 
                 U_SS=s$U_SS_list[[ind]], 
                 ellipses=TRUE,
                 gg.add=list(theme_bw()),
                 ...
        )
    }else if(s$clust_distrib=="skewNormal"){
        plot_DPMsn(z=s$data,
                   c=s$point_estim$c_est, 
                   i=ind+s$burnin, 
                   alpha=s$alpha[ind], 
                   U_SS=s$U_SS_list[[ind]], 
                   ellipses=TRUE,
                   gg.add=list(theme_bw()),
                   nbsim_dens=nbsim_densities,
                   ...
        )
    }else if(s$clust_distrib=="skewT"){
        plot_DPMst(z=s$data,
                   c=s$point_estim$c_est, 
                   i=ind+s$burnin, 
                   alpha=s$alpha[ind], 
                   U_SS=s$U_SS_list[[ind]], 
                   ellipses=TRUE,
                   gg.add=list(theme_bw()),
                   nbsim_dens=nbsim_densities,
                   ...
        )
    }
    cat("DONE!\n")
    
    if(hm){
        cat("Plotting heatmap of similarity (may take a few min)...\n")
        pheatmap(s$point_estim$similarity, scale="none", border_color=NA,
                 color=colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", 
                                          #"#9ECAE1", "#FEB24C", 
                                          "#FD8D3C", "#BD0026", "#800026"))(200), 
                 show_rownames=FALSE, show_colnames=FALSE, 
                 cluster_rows=TRUE, cluster_cols =TRUE, 
                 main="Posterior similarity matrix")
        cat("DONE! Wait for plot rendering...\n")
    }
    
}
