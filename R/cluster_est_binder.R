#'Gets a point estimate of the partition using Binder loss function
#'
#'@param c a list of vector of length \code{n}. \code{c[[j]][i]} is 
#'the cluster allocation of observation \code{i=1...n} at iteration 
#'\code{j=1...N}.
#'
#'@param thin integer indicating the thinning of the MCMC output. 
#'Default is \code{1}.
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
#'
#'@author Francois Caron, Boris Hejblum
#'
#'@export cluster_est_binder
#'
#'@references F. Caron, Y.W. Teh, T.B. Murphy. Bayesian nonparametric Plackett-Luce 
#' models for the analysis of preferences for college degree programmes. 
#' To appear in Annals of Applied Statistics, 2014.
#' http://arxiv.org/abs/1211.5037
#' 
#' D. B. Dahl. Model-Based Clustering for Expression Data via a 
#' Dirichlet Process Mixture Model, in Bayesian Inference for 
#' Gene Expression and Proteomics, K.-A. Do, P. Muller, M. Vannucci 
#' (Eds.), Cambridge University Press, 2006.
#'
#'@seealso \link{similarityMat}
#'

cluster_est_binder <- function(c, thin=1){
    n <- length(c[[1]])
    
    if(thin>1){
        select <- c(TRUE, rep(FALSE, thin-1))
        c <- c[select]
    }
    
    N <- length(c)

    
    #Non vectorized piece of code for reference
#     cost <- numeric(N)
#     similarity <- matrix(ncol=n, nrow=n)
#     for(i in 1:(n-1)){
#          for(j in (i+1):n){
#              similarity[i,j] <- 1/N*sum(unlist(lapply(c, "[", i)) == unlist(lapply(c, "[", j)))    
#              for(k in 1:N){
#                  cost[k] = cost[k] + abs(as.numeric(c[[k]][i]==c[[k]][j]) - similarity[i,j])
#              }
#         }
#     }
#     cost <- 2*cost
#     which.min(cost)
    
    
    vclust2mcoclust <- function(v){
        m <- sapply(v, FUN=function(x){x==v})
        return(m)
    }
    
    list_mcoclust <- lapply(c, vclust2mcoclust)
    
    similarityFast <- Reduce('+', list_mcoclust)/N
    list_cost <- lapply(list_mcoclust, function(m){abs(m-similarityFast)})
    costFast <- unlist(lapply(list_cost, sum))
    opt_ind <- which(costFast==min(costFast))
    opt_ind <- opt_ind[length(opt_ind)]
    c_est <- c[[opt_ind]]
    
    return(list("c_est"=c_est, "cost"=costFast, "similarity"=similarityFast, "opt_ind"=opt_ind))
}