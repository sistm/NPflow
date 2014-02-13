#' @keywords internal

plot_DPM <- function(z, U_mu, U_Sigma, m, c, i){
    fullCl <- which(m!=0)
    U_mu2plot <- U_mu[, fullCl]    
    zClusters <- factor(c, levels=as.character(fullCl), ordered=TRUE)
    #levels(zClusters) <- as.character(1:length(levels(zClusters)))
    
    z2plot <- cbind.data.frame("X"=z[1,],"Y"=z[2,],"Cluster"=zClusters)
    if(is.null(dim(U_mu2plot))){
        U2plot <- cbind.data.frame("X"=U_mu2plot[1],
                                   "Y"=U_mu2plot[2],
                                   "Cluster"=factor(as.character(fullCl), 
                                                    levels=as.character(fullCl), 
                                                    ordered=TRUE)
        )
        
    } else {
        U2plot <- cbind.data.frame("X"=U_mu2plot[1,],
                                   "Y"=U_mu2plot[2,],
                                   "Cluster"=factor(as.character(fullCl), 
                                                    levels=as.character(fullCl), 
                                                    ordered=TRUE)
        )     
    }
    
    p <- (ggplot(z2plot) 
          + geom_point(aes(x=X, y=Y, colour=Cluster, order=Cluster), data=z2plot)
          + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster), data=U2plot, shape=22, size=5)
          + ggtitle(paste("Iteration", i))
    )
    print(p)
}