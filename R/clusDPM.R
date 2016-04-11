#'Unsupervised clustering though Dirichlet Process Mixture models
#'
#'Wrapping function for mcmc sampling and partition estimation.
#'
#'@param data
#'
#'@param distrib
#'
#'@param ncores
#'
#'@seealso \code{\link{DPMpost}}
#'@seealso \code{\link{cluster_est_binder}}
#'
#'@export
#'
clusDPM <- function(data, distrib=c("gaussian", "skewnorm", "skewt"), ncores=1
){
    res <- list(data)
    class(res) <- "NPlow"

    #check hyperG0
}