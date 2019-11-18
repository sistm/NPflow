#'Scatterplot of flow cytometry data
#'
#'@param cytomatrix a \code{p x n} data matrix, of \code{n} cell observations measured over \code{p}
#'markers.
#'
#'@param dims2plot a vector of length at least 2, indicating of the dimensions to be plotted.
#'Default is \code{c(1, 2)}.
#'
#'@param gating an optional vector of length \code{n} indicating a known gating of the cells to be 
#'displayed. Default is \code{NULL} in which case no gating is displayed.
#'
#'@param scale_log a logical Flag indicating whether the data should be plotted on the log scale.
#'Default is \code{FALSE}.
#'
#'@param xlim a vector of length 2 to specify the x-axis limits. Only used if \code{dims2plot} is
#'of length 2Default is the data range.
#'
#'@param ylim a vector of length 2 to specify the y-axis limits. Only used if \code{dims2plot} is
#'of length 2. Default is the data range.
#'
#'@param gg.add
#'A list of instructions to add to the \code{ggplot2} instruction (see \code{\link[ggplot2]{gg-add}}).
#'Default is \code{list(theme())}, which adds nothing.
#'to the plot.
#'
#'@examples
#'
#' rm(list=ls())
#' #Number of data
#' n <- 500
#' #n <- 2000
#' set.seed(1234)
#' #set.seed(123)
#' #set.seed(4321)
#'
#' # Sample data
#' m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#'
#' sdev <- array(dim=c(2,2,4))
#' sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#' c <- rep(0,n)
#' z <- matrix(0, nrow=2, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' cytoScatter(z)
#'
#'@export
cytoScatter <- function(cytomatrix, dims2plot=c(1,2), gating=NULL,
                        scale_log=FALSE, xlim=NULL, ylim=NULL, gg.add=list(theme())){
  
  n <- ncol(cytomatrix)
  stopifnot(length(dims2plot)>1)
  stopifnot(is.matrix(cytomatrix))
  
  
  if(length(dims2plot)!=2){
    cytomatrix <- cytomatrix[dims2plot,]
    
    zDplot <- melt(cbind.data.frame("ID"=as.character(1:n), t(cytomatrix)),
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
    
    if(!is.null(gating)){
      stopifnot(length(gating)==n)
      zDplotfull <- cbind.data.frame(zDplotfull, "Gate"=rep(as.character(gating), nrow(zDplotfull)/length(gating)))
    }
    
    if(is.null(gating)){
      p <- (ggplot(zDplotfull, aes_string(x="X", y="Y"))
            + facet_grid(dimensionY~dimensionX, scales="free")
            + geom_point(colour="blue", data=zDplotfull, alpha=1, size=2/(0.3*log(n)))
            + stat_density2d(aes_string(fill = "..level.."), alpha=0.8, geom="polygon")
            + scale_fill_gradientn(colours=c("blue","green", "yellow", "red"),
                                   name="Density")
            + theme_bw()
            + ggtitle(paste(n, " cells", sep=""))
      )
    }else{
      p <- (ggplot(zDplotfull, aes_string(x="X", y="Y", group="Gate"))
            + facet_grid(dimensionY~dimensionX, scales="free")
            #+ stat_density2d(aes_string(fill = "..level.."), alpha=0.1, geom="polygon")
            + geom_point(aes_string(colour="Gate"), data=zDplotfull, alpha=0.5/(0.3*log(n)), size=2/(0.3*log(n)), shape=16)      
            + theme_bw()
            + ggtitle(paste(n, " cells", sep=""))
            + guides(colour=guide_legend(override.aes = list(alpha=1)))
      )
    }
    
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
    if(!is.null(gating)){
      stopifnot(length(gating)==n)
      data2plot <- cbind.data.frame(data2plot, "Gate"=as.character(gating))
    }
    
    if(is.null(colnames(data2plot)[1]) | is.null(colnames(data2plot)[1])){
      colnames(data2plot)[1:2] <- c("X", "Y")
    }
    xname <- colnames(data2plot)[1]
    yname <- colnames(data2plot)[2]
    p <- ggplot(data2plot, aes_string(x=xname, y=yname))

    if(is.null(gating)){
    p <- (p + geom_point(size=2/(0.3*log(n)), colour="blue")
          + stat_density2d(aes_string(fill = "..level.."), alpha=0.8, geom="polygon")
          + scale_fill_gradientn(colours=c("blue","green", "yellow", "red"),
                                 name="Density")
          + theme_bw()
          + ggtitle(paste(n, " cells", sep=""))
    )
    }else{
      p <- (p #+ stat_density2d(aes_string(fill = "..level.."), alpha=0.1, geom="polygon")
            + geom_point(aes_string(colour="Gate"), size=2/(0.3*log(n)), alpha=0.4)
            + theme_bw()
            + guides(colour=guide_legend(override.aes = list(alpha=1)))
            + ggtitle(paste(n, " cells", sep=""))
      )
    }
    
    
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