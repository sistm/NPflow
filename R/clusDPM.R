#'Unsupervised clustering though Dirichlet Process Mixture models
#'
#'Wrapping function for mcmc sampling and partition estimation.
#'
#'@param data a data.frame or matrix
#'
#'@param distribution a character string, either \code{"Normal"}, \code{"skewNormal"}, \code{"skewStudent"}.
#'
#'@seealso DPM_GibbsSampler_SkewN
#'@seealso cluster_est_binder
#'
#'@export
#'
clusDPM <- function(data,
                      distibution=c("Normal", "skewNormal", "skewStudent")
){
    if(distribution=="skewStudent"){
        stop("skew-Student distribution is not implemented yet\n")
    }
    
    TODO
  
    res <- list(data)
    class(res) <- "NPlow"
}