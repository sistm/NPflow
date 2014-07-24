#'Scatterplot of flow cytometry data
#'
#'@param cytomatrix
#'
#'@param dims2plot
#'
#'@param scale_log
#'
#'@param xlim
#'
#'@param ylim
#'
#'@export
#'
#'

cytoScatter <- function(cytomatrix, dims2plot=c(1,2), 
                      scale_log=FALSE, xlim=NULL, ylim=NULL, gg.add=list(theme())){
    if(length(dims2plot)!=2){
        
        n <- ncol(cytomatrix)
        cytomatrix <- cytomatrix[dims2plot,]
        
        zDplot <- melt(cbind.data.frame("ID"=as.character(1:n), 
                                        t(cytomatrix)
        ),
        id.vars=c("ID"), 
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
        
        p <- (ggplot(zDplotfull, aes(x=X, y=Y)) 
              + facet_grid(dimensionY~dimensionX, scales="free")
              + geom_point(colour="blue",
                           data=zDplotfull, alpha=1, size=2/(0.3*log(n)))
              + stat_density2d(aes(fill = ..level..), alpha=0.8, geom="polygon")  
              + scale_fill_gradientn(colours=c("blue","green", "yellow", "red"), 
                                     name="Density")
              + theme_bw()
              + ggtitle(paste(n, " cells",
                              sep=""))
        )
        
    }else{
        
        if(is.character(dims2plot)){
            if(length(which(dims2plot%in%rownames(cytomatrix)))!=2){
                stop("'dims2plot' not in rownames of 'cytomatrix'")
            }
            dims2plot[1] <- which(rownames(cytomatrix)==dims2plot[1])
            dims2plot[2] <- which(rownames(cytomatrix)==dims2plot[2])
            dims2plot <- as.numeric(dims2plot)
        }
        
        data2plot <- data.frame(t(cytomatrix[dims2plot,]))
        if(!is.null(colnames(data2plot)[dims2plot[1]]) && 
               !is.null(colnames(data2plot)[dims2plot[2]])){
            xname <- colnames(data2plot)[dims2plot[1]]
            yname <- colnames(data2plot)[dims2plot[2]]
            p <- ggplot(data2plot, aes_string(x=xname, y=yname))
        }else{
            colnames(data2plot) <- c("X", "Y")
            p <- ggplot(data2plot, aes_string(x=X, y=Y))
        }
        
        p <- (p + geom_point(size=0.4, colour="blue")  
              + stat_density2d(aes(fill = ..level..), alpha=0.8, geom="polygon")  
              + scale_fill_gradientn(colours=c("blue","green", "yellow", "red"), 
                                     name="Density")
              + theme_bw()
        )
        
        
        if (!is.null(xlim)){
            if(scale_log){
                p <- (p + scale_x_log10(limits=xlim, colnames(data2plot)[dims2plot[1]]))
                if (!is.null(ylim)){
                    p <- (p + scale_y_log10(limits=ylim, colnames(data2plot)[dims2plot[2]]))
                }
            }
            else{
                p <- (p + xlim(xlim))
                if (!is.null(ylim)){
                    p <- (p + ylim(ylim))
                }
            }
        }
        if (!is.null(ylim) && is.null(xlim)){
            if(scale_log){
                p <- (p + scale_y_log10(limits=ylim, colnames(data2plot)[dims2plot[2]])
                      + scale_x_log10(colnames(data2plot)[dims2plot[1]])
                )
            }
            else{
                p <- (p + ylim(ylim))
            }
        }
        
        if(scale_log && is.null(ylim) && is.null(xlim)){
            p <- (p + scale_x_log10(colnames(data2plot)[dims2plot[1]])
                  + scale_y_log10(colnames(data2plot)[dims2plot[2]])
            )
        }
    }
    
    for (a in gg.add) {
        p <- p + a
    }
    print(p)
}