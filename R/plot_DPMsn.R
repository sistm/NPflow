#' @author Boris Hejblum
#'
#' @export plot_DPMsn

plot_DPMsn <- function(z, U_xi, U_psi, U_Sigma, m, c, i, alpha="?", U_SS=NULL,
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
    
    if(is.list(U_xi)){
        U_mu2plot <- matrix(0, nrow=p, ncol=length(fullCl))
        U_Sigma2plot <- array(0, dim=c(p, p, length(fullCl)))
        for(i in 1:length(fullCl)){
            k <- as.character(fullCl[i])
            U_mu2plot[,i] <- U_xi[[k]][dims2plot]+U_psi[[k]][dims2plot]*mean_sn01
            colnames(U_mu2plot) <- fullCl
            rownames(U_mu2plot) <- rownames(z)
            U_Sigma2plot[, , i] <- U_Sigma[[k]][dims2plot, dims2plot]
        }
    }else{
        U_mu2plot <- U_xi[, fullCl, drop=FALSE] + U_psi[, fullCl, drop=FALSE]*mean_sn01
        rownames(U_mu2plot) <- rownames(z)
        U_Sigma2plot <- U_Sigma[, , fullCl]
        U_xi2plot <- U_xi[, fullCl, drop=FALSE]
        rownames(U_xi2plot) <- rownames(z)
    }
    U_SS2plot <- U_SS#[fullCl]
    zClusters <- factor(c, levels=as.character(fullCl), ordered=TRUE)
    #levels(zClusters) <- as.character(1:length(levels(zClusters)))
    
    expK <- ifelse(is.numeric(alpha), round(alpha*(digamma(alpha+n)-digamma(alpha))), NA)
    alpha2print <- ifelse(is.numeric(alpha), formatC(alpha, digits=2), alpha)
    
    library(ellipse)
    
    if(p>2){
        zDplot <- melt(cbind.data.frame("ID"=as.character(1:n), 
                                        t(z),
                                        "Cluster"=zClusters
        ),
        id.vars=c("ID", "Cluster"), 
        variable.name = "dimensionX",
        value.name="X"
        )
        zDplotfull <- zDplot 
        zDplotfull$Y <- zDplot$X
        zDplotfull$dimensionY <- zDplot$dimensionX
        
        lev <- as.character(1:length(levels(zDplot$dimensionX)))
        for(l in 2:length(lev)){
            move <- which(as.numeric(zDplot$dimensionX)<l)
            zDplottemp <- rbind.data.frame(zDplot[-move,], zDplot[move,])
            zDplottemp$Y <- zDplot$X
            zDplottemp$dimensionY <- zDplot$dimensionX
            zDplotfull <- rbind.data.frame(
                zDplotfull, zDplottemp)
        }
        
        UDplot <- melt(cbind.data.frame(t(U_mu2plot),
                                        "Cluster"=factor(as.character(fullCl), 
                                                         levels=as.character(fullCl), 
                                                         ordered=TRUE)
        ),
        id.vars=c("Cluster"), 
        variable.name = "dimensionX",
        value.name="X"
        )
        UDplotfull <- UDplot 
        UDplotfull$Y <- UDplotfull$X
        UDplotfull$dimensionY <- UDplotfull$dimensionX
        
        lev <- levels(UDplotfull$dimensionX)
        for(l in 2:length(lev)){
            move <- which(as.numeric(UDplotfull$dimensionX)<l)
            UDplottemp <- rbind.data.frame(UDplotfull[-move,], UDplotfull[move,])
            UDplottemp$Y <- UDplotfull$X
            UDplottemp$dimensionY <- UDplotfull$dimensionX
            UDplotfull <- rbind.data.frame(
                UDplotfull, UDplottemp)
        }
        
        #         ellipse95 <- data.frame()
        #         for(dx in 1:p){
        #             for(dy in 1:p){
        #                 for(g in 1:length(fullCl)){
        #                     glabel <- levels(zClusters)[g]
        #                     U_corr2plot_g <- cov2cor(U_Sigma2plot[c(dx,dy),c(dx,dy),g])
        #                     ellipse95 <- rbind(ellipse95, 
        #                                cbind(as.data.frame(ellipse(U_corr2plot_g, 
        #                                                            scale=sqrt(diag(U_Sigma2plot[c(dx,dy),c(dx,dy),g])), 
        #                                                            centre=U_mu2plot[c(dx,dy),g])),
        #                                      Cluster=as.character(glabel), 
        #                                      dimensionX=as.character(dx), 
        #                                      dimensionY=as.character(dy)
        #                                ))
        #                 }
        #             }
        #         }
        
        p <- (ggplot(zDplotfull) 
              + facet_grid(dimensionY~dimensionX, scales="free")
              + geom_point(aes(x=X, y=Y, colour=Cluster, order=Cluster), 
                           data=zDplotfull, alpha=1, size=2/(0.3*log(n)))
              #               + geom_polygon(aes(x=x, y=y, fill=Cluster, colour=Cluster, order=Cluster), 
              #                              data=ellipse95, size=0.5, linetype=2, colour="black", alpha=.3)
              + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster),
                           data=UDplotfull, shape=22, size=5/(0.3*log(n)))
              + ggtitle(paste(n, " obs.",
                              "\niteration ", i, " : ", 
                              length(fullCl)," clusters",
                              "\nexpected number of clusters: ", expK,
                              " (alpha = ", alpha2print, ")",
                              sep=""))
        )
    }else{
        z2plot <- cbind.data.frame("D1"=z[1,],"D2"=z[2,],"Cluster"=zClusters)
        U2plot <- cbind.data.frame("D1"=U_mu2plot[1,],
                                   "D2"=U_mu2plot[2,],
                                   "Cluster"=factor(as.character(fullCl), 
                                                    levels=as.character(fullCl), 
                                                    ordered=TRUE)
        )  
        xi2plot <- cbind.data.frame("D1"=U_xi2plot[1,],
                                    "D2"=U_xi2plot[2,],
                                    "Cluster"=factor(as.character(fullCl), 
                                                     levels=as.character(fullCl), 
                                                     ordered=TRUE)
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
        if(ellipses){
            simuDens <- NULL
            for(g in 1:length(fullCl)){
                glabel <- levels(zClusters)[g]
                gind <- as.numeric(glabel)
                ltnz <- rtruncnorm(5000, a=0)
                eps <- matrix(rnorm(2*5000), ncol=2)%*%chol(U_Sigma[,,gind])
                simuDenstemp <- data.frame("D1"=U_xi[1,gind]+U_psi[1,gind]*ltnz+eps[,1], 
                                           "D2"=U_xi[2,gind]+U_psi[2,gind]*ltnz+eps[,2], 
                                           "Cluster"=rep(glabel, 5000))
                simuDens <- rbind.data.frame(simuDens, simuDenstemp)
            }
            
            p <- (p 
                  + geom_density2d(data=simuDens, aes(x=D1,y=D2, colour=Cluster, linetype="1"))
                  + scale_linetype_manual(values=c(1),
                                          labels=c("simulations derived\n from sampled xi & psi"),
                                          name="Density contour", limits=c(1))
            )
        }
    }
    for (a in gg.add) {
        p <- p + a
    }
    print(p)
}