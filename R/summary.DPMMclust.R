#'Summarizing Dirichlet Process Mixture Models
#'
#'Summary methods for \code{DPMMclust} objects.
#'
#'@param object a \code{DPMMclust} object.
#'
#'@param burnin integer giving the number of MCMC iterations to burn (defaults is half)
#'
#'@param thin integer giving the spacing at which MCMC iterations are kept.
#'Default is \code{1}, i.e. no thining.
#'
#'@param lossFn character string specifying the loss function to be used.
#'Either "F-measure" or "Binder" (see Details). Default is "F-measure".
#'
#'@param gs optional vector of length \code{n} containing the gold standard
#'partition of the \code{n} observations to compare to the point estimate
#'
#'@param posterior_approx logical flag whether a parametric approximation of the posterior should be
#'computed. Default is \code{FALSE}
#'
#'@param ... further arguments passed to or from other methods
#'
#'@return a \code{list} containing the following elements:
#'  \itemize{
#'      \item{\code{nb_mcmcit}:}{ an integer giving the value of \code{m}, the number of retained 
#'      sampled partitions, i.e. \code{(N - burnin)/thin}}
#'      \item{\code{burnin}:}{ an integer passing along the \code{burnin} argument}
#'      \item{\code{thin}:}{ an integer passing along the \code{thin} argument}
#'      \item{\code{lossFn}:}{ a character string passing along the \code{lossFn} argument}
#'      \item{\code{clust_distrib}:}{ a character string passing along the \code{clust_distrib} argument }
#'      \item{\code{point_estim}:}{ a \code{list} containing: \itemize{
#'          \item{\code{c_est}:}{ a vector of length \code{n}containing the point estimated clustering for each observations}
#'          \item{\code{cost}:}{ a vector of length \code{m} containing the cost of each sampled partition}
#'          \item{\code{Fmeas}:}{ if \code{lossFn} is \code{'F-measure'}, the \code{m x m} matrix of total F-measures for each pair of sampled partitions}
#'          \item{\code{opt_ind}:}{ the index of the point estimate partition among the \code{m} sampled}
#'      }}
#'      \item{\code{loss}:}{ the loss for the point estimate. \code{NA} if \code{lossFn} is not \code{'Binder'}}
#'      \item{\code{param_posterior}:}{ a list containing the parametric approximation of the posterior,
#'      suitable to be plugged in as prior for a new MCMC algorithm run}
#'      \item{\code{mcmc_partitions}:}{ a list containing the \code{m} sampled partitions}
#'      \item{\code{alpha}:}{ a vector of length \code{m} with the values of the \code{alpha} DP parameter}
#'      \item{\code{index_estim}:}{ the index of the point estimate partition among the \code{m} sampled}
#'      \item{\code{hyperG0}:}{ a list passing along the prior, i.e. the \code{hyperG0} argument}
#'      \item{\code{logposterior_list}:}{ a list of length \code{m} containing the logposterior and its decomposition, for each sampled partition}
#'      \item{\code{U_SS_list}:}{ a list of length \code{m} containing the containing the lists of sufficient statistics for all the mixture components,
#'      for each sampled partition}
#'      \item{\code{data}: a \code{d x n} matrix containing the clustered data}
#'  }
#'  
#'
#'@details The cost of a point estimate partition is calculated using either a pairwise
#' coincidence loss function (Binder), or 1-Fmeasure (F-measure).
#' 
#'The number of retained sampled partitions is \code{m = (N - burnin)/thin}
#'
#'@author Boris Hejblum
#'
#'@export
#'
#'@importFrom pheatmap pheatmap
#'
#'@importFrom fastcluster hclust
#'
#'@seealso \code{\link{similarityMat}} \code{\link{similarityMatC}}
#'

summary.DPMMclust <- function(object, burnin=0, thin=1, gs=NULL, lossFn="F-measure",
                              posterior_approx=FALSE, ...){
  
  x_invar <- burn.DPMMclust(object, burnin = burnin, thin=thin)
  
  
  
  #if(!posterior_approx){
  if(lossFn == "F-measure"){
    point_estim <- cluster_est_Fmeasure(x_invar$mcmc_partitions,
                                        logposterior = sapply(x_invar$logposterior_list, sum))
  }else if(lossFn == "Binder"){
    point_estim <- cluster_est_binder(x_invar$mcmc_partitions,
                                      logposterior = sapply(x_invar$logposterior_list, sum))
  }else if(lossFn == "MBinderN"){
    if (x_invar$clust_distrib!="gaussian"){stop("Clusters distribution must be Gaussian")}
    point_estim <- cluster_est_Mbinder_norm(x_invar$mcmc_partitions,x_invar$listU_mu,
                                            x_invar$listU_Sigma,
                                            logposterior = sapply(x_invar$logposterior_list, sum), ...)
  }else{
    stop("Specified loss function not available.\n
         Specify either 'F-measure' or 'Binder' for the lossFn argument.")
  }
  
  index_estim <- point_estim$opt_ind
  loss <- NA
  if(!is.null(gs)){
    loss <- evalClustLoss(c=point_estim$c_est, gs=gs, lossFn=lossFn)
  }
  #}
  #point_estim=NULL
  #index_estim <- NA
  #loss <- NA
  #Posterior approximation
  K<-length(unique(point_estim$c_est))
  if(posterior_approx){
    if(K>1){
      param_post <- postProcess.DPMMclust(x_invar, plot=FALSE, K=K, ...)
    }else{
      param_post <- postProcess.DPMMclust(x_invar, K=1, ...)
    }
  }
  else{
    param_post <- NULL
  }
  
  s <- c(x_invar, list("burnin"=burnin,
                       "thin"=thin,
                       "point_estim"=point_estim,
                       "loss"=loss,
                       "index_estim"=index_estim,
                       "param_posterior"=param_post))
  class(s) <- "summaryDPMMclust"
  
  invisible(s)
}


#' Methods for a summary of a \code{DPMMclust} object
#'@rdname methods.summaryDPMMclust
#'@aliases summaryDPMMclust print.summaryDPMMclust plot.summaryDPMMclust
#'@param x a \code{summaryDPMMclust} object.
#'@param hm logical flag to plot the heatmap of the similarity matrix.
#'Default is \code{FALSE}.
#'@param nbsim_densities the number of simulated observations to be used
#'to plot the density lines of the clusters in the point estimate partition plot
#'@param ... further arguments passed to or from other methods
#'@author Boris Hejblum
#'@export
print.summaryDPMMclust <- function(x,...){
  
  cat(class(x), "object with", x$nb_mcmcit, "observations sampled from the posterior:\n")
  cat(rep("-", 72), "\n", sep="")
  #cat(rep("-",40),"\n", sep="")
  cat("Burnin:", x$burnin, "MCMC iterations discarded\n\n")
  cat("Point estimate of the partition with a ", x$clust_distrib," mixture:\n", sep="")
  t<-table(x$point_estim$c_est)
  t2print <- paste(formatC(as.vector(t/length(x$point_estim$c_est)),
                           digits=2),"%", sep="")
  names(t2print) <- names(t)
  print(t2print, quote=F)
  if(!is.na(x$loss)){
    cat("\nLoss of the point estimate partition compared to the reference clustering provided:", x$loss, "\n")
    cat(rep("-",40),"\n", sep="")
  }
}


#'@rdname methods.summaryDPMMclust
#'@param gg.add a list of instructions to add to the \code{ggplot2} instruction (see 
#'\code{\link[ggplot2]{gg-add}}). Default is \code{list(theme())}, which adds nothing to the plot.
#'@param hm_subsample a integer designating the number of observations to use when plotting the heatmap. 
#'Used only if \code{hm} is \code{TRUE}. #'Default is \code{NULL} in which no subsampling is done and 
#'all observations are plotted.
#'@param hm_order_by_clust logical flag indicating whether observations should be ordered according to
#'the point estimate first. Used only if \code{hm} is \code{TRUE}. Default is \code{TRUE}.
#'@importFrom stats dist
#'@importFrom grDevices colorRampPalette
#'@export
plot.summaryDPMMclust <- function(x, hm=FALSE, nbsim_densities=5000, 
                                  hm_subsample=NULL, hm_order_by_clust=TRUE, 
                                  gg.add=list(theme_bw()),...){
  
  if(length(x$logposterior_list[[1]])>1){
    plot_ConvDPM(x, shift=x$burnin, thin=x$thin)
  }
  
  ind <- x$index_estim
  
  cat("Plotting point estimate (may take a few sec)... ")
  if(x$clust_distrib=="gaussian"){
    plot_DPM(z=x$data,
             c=x$point_estim$c_est,
             i=ind+x$burnin,
             alpha=x$alpha[ind],
             U_SS=x$U_SS_list[[ind]],
             ellipses=TRUE,
             gg.add=gg.add,
             ...
    )
  }else if(x$clust_distrib=="skewnorm"){
    plot_DPMsn(z=x$data,
               c=x$point_estim$c_est,
               i=ind+x$burnin,
               alpha=x$alpha[ind],
               U_SS=x$U_SS_list[[ind]],
               ellipses=TRUE,
               gg.add= gg.add,
               nbsim_dens=nbsim_densities,
               ...
    )
  }else if(x$clust_distrib=="skewt"){
    plot_DPMst(z=x$data,
               c=x$point_estim$c_est,
               i=ind+x$burnin,
               alpha=x$alpha[ind],
               U_SS=x$U_SS_list[[ind]],
               ellipses=TRUE,
               gg.add=gg.add,
               nbsim_dens=nbsim_densities,
               nice=TRUE,
               ...
    )
  }
  cat("DONE!\n")
  
  if(hm){
    
    if(is.null(x$point_estim$similarity)){
      warning("\nThe similarity matrix makes more sense when using the 'Binder' loss function...")
      cat("Estimating posterior similarity matrix (this may take some time, complexity in O(n^2))...")
      x$point_estim$similarity <- similarityMat_nocostC(do.call(cbind, x$mcmc_partitions))$similarity
      cat("Done!\n")
    }
    
    cat("Plotting heatmap of similarity (may take a few min):\n")
    if(!is.null(hm_subsample)){
      if(length(hm_subsample)>1){
        select <- hm_subsample
      }else{
        select <- sample(1:ncol(x$point_estim$similarity), size=hm_subsample)
      }
      x$point_estim$similarity <- x$point_estim$similarity[select, select]
      x$point_estim$c_est <- x$point_estim$c_est[select]
    }
    
    
    if(hm_order_by_clust){
      cat(" computing pairwise distances...")
      prop <- table(x$point_estim$c_est)
      clusters_ordered <- names(prop)[order(prop, decreasing=TRUE)]
      
      ord_index <- list()
      for(k in clusters_ordered){
        index <- which(x$point_estim$c_est==as.integer(k))
        if(length(index)>1){
          dist_mat <- stats::dist(x$point_estim$similarity[index, index], method = "euclidean")
          tree <- fastcluster::hclust(dist_mat, method = "complete")
          ord_index[[k]] <- index[tree$order]
        }else{
          ord_index[[k]] <- index
        }
      }
      cat("DONE!\n ordering samples now...")
      ord_final <- do.call(c, ord_index)
      cat("DONE!\n")
    }else{
      cat(" computing pairwise distances...")
      dist_mat <- stats::dist(x$point_estim$similarity, method = "euclidean")
      cat("DONE!\n ordering samples now...")
      tree <- fastcluster::hclust(dist_mat, method = "complete")
      ord_final <- tree$order
      cat("DONE!\n")
    }
    
    my_annot_row <- cbind.data.frame("Clustering from point estimate"=factor(x$point_estim$c_est))
    rownames(x$point_estim$similarity) <- as.character(1:nrow(x$point_estim$similarity)) 
    colnames(x$point_estim$similarity) <- as.character(1:ncol(x$point_estim$similarity)) 
    
    pheatmap::pheatmap(x$point_estim$similarity[ord_final, ord_final], scale="none", border_color=NA,
                       color=grDevices::colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF",
                                                           #"#9ECAE1", "#FEB24C",
                                                           "#FD8D3C", "#BD0026", "#800026"))(200),
                       show_rownames=FALSE, show_colnames=FALSE,
                       cluster_rows=FALSE, cluster_cols = FALSE,
                       #annotation_row = my_annot_row[ord_index, , drop=FALSE],
                       annotation_col = my_annot_row[ord_final, , drop=FALSE],
                       annotation_names_col = FALSE,
                       main="Posterior similarity matrix\n")
    cat("Almost over, wait for plot rendering now...\n\n")
  }
}
