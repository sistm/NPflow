#'Unsupervised clustering though Dirichlet Process Mixture models
#'
#'Wrapping function for mcmc sampling and partition estimation.
#'
#'@param data
#'
#'@param distrib
#'
#'@seealso \code{\link{DPMpost}}
#'@seealso \code{\link{cluster_est_binder}}
#'
#'@export
#'
clusDPM <- function(data, distib=c("gaussian", "skewnorm", "skewt")
){
    res <- list(data)
    class(res) <- "NPlow"

    #check hyperG0
}