#'Unsupervised clustering though Dirichlet Process Mixture models
#'
#'Wrapping function for mcmc sampling and partition estimation.
#'
#'@param data
#'
#'@param distribution
#'
#'@seealso DPM_GibbsSampler_SkewN
#'@seealso cluster_est_binder
#'
#'@export clusDPM
#'
clusDPM <- function(data,
                      distibution=c("Normal", "skewNormal", "skewStudent")
){
    if(distribution=="skewStudent"){
        stop("Skew-Student distribution is not implemented yet\n")
    }
    
    res <- list(data, )
    
    class(res) <- "NPlow"
}