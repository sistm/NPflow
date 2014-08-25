#'Compute F-measure
#'
#'@name Fmeasure
#'
#'@param pred
#'
#'@param ref
#'
#'@param select_ref
#'
#'@return a list: 
#'  \itemize{
#'      \item{\code{Fsum}:}{}
#'      \item{\code{Fmeasures}:}{}
#'  }
#'
#'@author Boris P. Hejblum
#'
#'@references Aghaeepour, N., Finak, G., Hoos, H., Mosmann, 
#'T. R., Brinkman, R., Gottardo, R., & Scheuermann, R. H. 
#'Critical assessment of automated flow cytometry data analysis techniques. 
#'Nature Methods, 10(3): 228-38, 2013.
#'
#'@export Fmeasure
#'
#'
#'
Fmeasure <- function(pred, ref, select_ref=NULL){
    #@example 
    #ftest <- Fmeasure(pred=s$point_estim$c_est, ref=reference.clusters[[1]], select_ref =order(reference.scores[[1]], decreasing=TRUE)[1:35])
    K <- unique(pred)
    K <- K[order(K)]
    C <- unique(ref)
    C <- C[order(C)]
    if(!is.null(select_ref)){
        C <- select_ref
    }
    m <- length(K)
    n <- length(C)
    
    M <- matrix(NA, nrow=n, ncol=m)
    Pr <- matrix(NA, nrow=n, ncol=m)
    Re <- matrix(NA, nrow=n, ncol=m)
    Fmat <- matrix(NA, nrow=n, ncol=m)
    C_card <- rep(NA, n)
    K_card <- rep(NA, m)
    for(i in 1:n){
        C_card[i] <- length(which(ref==C[i]))
    }
    for(j in 1:m){
        K_card[j] <- length(which(pred==K[j]))
    }
    
    
    for(i in 1:n){
        C_card[i] <- length(which(ref==C[i]))
       for(j in  1:m){
           M[i,j] <- length(which(ref==C[i] & pred==K[j]))
           Pr[i,j] <- M[i,j]/K_card[j]
           Re[i,j] <- M[i,j]/C_card[i] #length(which(c==C[i]))
           if(Pr[i,j]+Re[i,j]==0){
               Fmat[i,j] <- 0
           }else{
               Fmat[i,j] <-2*Pr[i,j]*Re[i,j]/(Pr[i,j]+Re[i,j])
           }
       }
    } 
    
    Ffinal <- apply(X=Fmat, MARGIN=1, max)
    names(Ffinal) <- C
    Fsum <- sum(Ffinal*C_card/sum(C_card))
    
    return(list("Fsum"=Fsum, "Fmeasures"=Ffinal)) 
}