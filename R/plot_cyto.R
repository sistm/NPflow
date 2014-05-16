plot_cyto <- function(cytomatrix, dims2plot=c(1,2)){
    if(length(dims2plot)!=2){
        stop("'dims2plot' must be of length 2")
    }
    
    if(is.character(dims2plot)){
        if(length(which(colnames(dims2plot%in%cytomatrix)))!=2){
            stop("'dims2plot' not in colnames of 'cytomatrix'")
        }
        dims2plot[1] <- which(colnames(cytomatrix)==dims2plot[1])
        dims2plot[2] <- which(colnames(cytomatrix)==dims2plot[2])
    }
    
    data2plot <- data.frame(cytomatrix[,c(dims2plot[1],dims2plot[2])])
    if(!is.null(colnames(data)[dims2plot[1]]) && !is.null(colnames(data)[dims2plot[2]])){
        xname <- colnames(data)[dims2plot[1]]
        yname <- colnames(data)[dims2plot[2]]
        p <- ggplot(data2plot, aes_string(x=xname, y=yname))
    }else{
        colnames(data2plot) <- c("X", "Y")
        p <- ggplot(data2plot, aes_string(x=X, y=Y))
    }
    
    (p + geom_point(size=0.4, colour="blue")  
     + stat_density2d(aes(fill = ..level..), alpha=0.8, geom="polygon")  
     + scale_x_log10(limits=c(25,50000), "CD27")
     + scale_y_log10(limits=c(30,20000), "CD20")
     + scale_fill_gradientn(colours=c("blue","green", "yellow", "red"), 
                            name="Density")
     + theme_bw()
    )
}