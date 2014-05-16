#' @author Boris Hejblum
#'
#' @export plot_DPMsn_test1K

plot_DPMsn_test1K <- function(z, U_xi, U_psi, U_Sigma, m, c, i, alpha="?", U_SS=NULL,
                       dims2plot=1:nrow(z),
                       ellipses=ifelse(length(dims2plot)<3,TRUE,FALSE),
                       gg.add=list(theme())){ 
    library(ellipse)
    library(reshape2)
    mean_sn01 <- (dnorm(0)-dnorm(Inf))/(pnorm(Inf)-pnorm(0))
    
    z <- z[dims2plot,]
    
    
    n <- ncol(z)
    p <- nrow(z)
    
    fullCl <- which(m!=0)
    

        U_mu2plot <- U_xi + U_psi*mean_sn01
        rownames(U_mu2plot) <- rownames(z)
        U_Sigma2plot <- U_Sigma
        U_xi2plot <- U_xi
        rownames(U_xi2plot) <- rownames(z)
    
    U_SS2plot <- U_SS#[fullCl]
    zClusters <- 1
    #levels(zClusters) <- as.character(1:length(levels(zClusters)))
    
    expK <- ifelse(is.numeric(alpha), round(sum(alpha/(alpha+1:(n-1)))), NA)
    alpha2print <- ifelse(is.numeric(alpha), formatC(alpha, digits=2), alpha)
    
    
        z2plot <- cbind.data.frame("D1"=z[1,],"D2"=z[2,],"Cluster"=as.factor(zClusters))
    
            U2plot <- cbind.data.frame("D1"=U_mu2plot[1],
                                       "D2"=U_mu2plot[2],
                                       "Cluster"="1"
            )  
            xi2plot <- cbind.data.frame("D1"=U_xi2plot[1],
                                        "D2"=U_xi2plot[2],
                                        "Cluster"="1"
            )   
    
        
        
        p <- (ggplot(z2plot) 
              
              + geom_point(aes(x=D1, y=D2, colour=Cluster, order=Cluster), 
                           data=z2plot)
              + geom_point(aes(x=D1, y=D2, fill=Cluster, order=Cluster, shape="22"),
                           data=U2plot, size=5)
              + geom_point(aes(x=D1, y=D2, fill=Cluster, order=Cluster, shape="23"),
                           data=xi2plot, size=5)
              + ggtitle(paste(n, " obs.",
                              "\niteration ", i, " : ", 
                              length(fullCl)," clusters",
                              "\nexpected number of clusters: ", expK,
                              " (alpha = ", alpha2print, ")",
                              sep=""))
              + scale_fill_discrete(guide=FALSE)
              + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6)))
              
        )
        
        #empirical mean of the clusters
        zmean2plot<- cbind.data.frame(D1=tapply(X=z2plot[,1], INDEX=z2plot$Cluster, FUN=mean),
                                      D2=tapply(X=z2plot[,2], INDEX=z2plot$Cluster, FUN=mean)
        )
        zmean2plot <- cbind.data.frame(zmean2plot, Cluster=rownames(zmean2plot))
        p <- (p + geom_point(aes(x=D1, y=D2, fill=Cluster, order=Cluster, shape="24"), 
                             data=zmean2plot, size=5)
              + scale_shape_manual(values=c(24,22,23), 
                                      labels=c("observed mean", "sampled mean", "xi param"), 
                                      name="", limits=c(24,22,23))
              )
    for (a in gg.add) {
        p <- p + a
    }    
    print(p)
}