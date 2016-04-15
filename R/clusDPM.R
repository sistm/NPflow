#'Unsupervised clustering though Dirichlet Process Mixture models
#'
#'Wrapping function for mcmc sampling and partition estimation.
#'
#'@param data a matrix data of dimension \code{n x p}.
#'
#'@param distrib the distribution used for the clustering. Current possibilities are
#'\code{"gaussian"}, \code{"skewnorm"} and \code{"skewt"}.
#'
#'@param ncores number of cores to use.
#'
#'@return an object of class \code{NPflow}. TODO
#'
#'@seealso \code{\link{DPMpost}}
#'@seealso \code{\link{cluster_est_binder}}
#'
#'@export
#'
clusDPM <- function(data, distrib=c("gaussian", "skewnorm", "skewt"), ncores=1){
  
  if(ncores < parallel::detectCores()){
    stop("Number of cores ")
  }
  
  res <- list(data)
  class(res) <- "NPlow"
  
  #check hyperG0
}