#'post-processing Dirichlet Process Mixture Models results to get 
#'a mixture distribution of the posterior locations
#'
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin integer giving the number of MCMC iterations to burn (defaults is half)
#'
#'@param thin integer giving the spacing at which MCMC iterations are kept. 
#'Default is \code{1}, i.e. no thining.
#'
#'@param lossFn character string specifying the loss function to be used.
#'Either "F-measure" or "Binder" (see Details). Default is "F-measure".
#'
#'@param gs optionnal vector of length \code{n} containing the gold standard 
#'partition of the \code{n} observations to compare to the point estimate
#'
#'@param ... further arguments passed to or from other methods
#'
#'@return a \code{list}: 
#'  \itemize{
#'      \item{\code{burnin}:}{an integer passing along the \code{burnin} argument}
#'      \item{\code{thin}:}{an integer passing along the \code{thin} argument}
#'      \item{\code{lossFn}:}{a character string passing along the \code{lossFn} argument}
#'      \item{\code{point_estim}:}{}
#'      \item{\code{loss}:}{}
#'      \item{\code{index_estim}:}{}
#'  }
#'
#'@details The cost of a point estimate partition is calculated using either a pairwise
#' coincidence loss function (Binder), or 1-Fmeasure (F-measure).
#'
#'@author Boris Hejblum
#'
#'@export 
#'
#'@importFrom gplots heatmap.2
#'
#'@seealso \link{similarityMat, summary.DPMMclust}
#'

postProcess.DPMMclust <- function(x, burnin=0, thin=1, gs=NULL, lossFn="F-measure", ...){
    
    x_invar <- burn.DPMMclust(x, burnin = burnin, thin=thin)
    
    S_final <- list()
    m_final <- list()
    for(i in 1:length(x_invar$U_SS_list)){
        xi_list <- sapply(x_invar$U_SS_list[[i]], "[", "xi")
        psi_list <- sapply(x_invar$U_SS_list[[i]], "[", "psi")
        m_final <- c(m_final, 
                     mapply(FUN=function(v1,v2){c(v1, v2)}, v1=xi_list, 
                            v2=psi_list, SIMPLIFY = FALSE)
        )
        
        S_list <- sapply(x_invar$U_SS_list[[i]], "[", "S")
        B_list <- sapply(x_invar$U_SS_list[[i]], "[", "B")
        S_final <- c(S_final, 
                                mapply(FUN=function(M1,M2){M1%x%M2}, M1=B_list, M2=S_list, SIMPLIFY = FALSE)
        )
        
    }
    
    
}
