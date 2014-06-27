#'Evaluate the loss of a point estimate of the partition compared to a gold standard through Binder loss function.
#'
#'The cost of a point estimate partition is calculated using a pairwise
#' coincidence loss function.
#'
#'@param c vector of length \code{n} containing the estimated partition 
#'of the  \code{n} observations.
#'
#'@param gs vector of length \code{n} containing the  gold standard 
#'partition of the  \code{n} observations.
#'
#'@param a penalty for wrong coclustering in \code{c} compared to code{gs}. Defaults is 1.
#'
#'@param b penalty for missed coclustering in \code{c} compared to code{gs}. Defaults is 1. 
#'
#'@return the cost of the point estimate \code{c} in regard of the 
#'gold standard \code{gs}.
#'
#'@author Boris Hejblum
#'
#'@export evalClust_Binder
#'
#'@references J.W. Lau & P.J. Green. Bayesian Model-Based Clustering 
#'Procedures, Journal of Computational and Graphical Statistics, 
#'16(3): 526–558, 2007.
#'
#' D. B. Dahl. Model-Based Clustering for Expression Data via a 
#' Dirichlet Process Mixture Model, in Bayesian Inference for 
#' Gene Expression and Proteomics, K.-A. Do, P. Muller, M. Vannucci 
#' (Eds.), Cambridge University Press, 2006.
#'
#'@seealso \link{similarityMat}, \link{clust_est_Binder}
#'

evalClust_Binder <- function(c, gs, a=1, b=1){
    n <- length(c)
    
    if(length(gs)!=n){
        stop("'c' and 'gs' arguments have not the same length")
    }
    
    c_coclust <- sapply(c, FUN=function(x){x==c})
    gs_coclust <- sapply(gs, FUN=function(x){x==gs})
    
    dif <- c_coclust-gs_coclust
    dif[which(dif=1)] <- b
    dif[which(dif=-1)] <- a
    
    loss <- sum(dif)
  
    
    return(loss)
}