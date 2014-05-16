\name{DPM_GibbsSampler_SkewN}
\alias{DPM_GibbsSampler_SkewN}
\title{Slice Sampling of Dirichlet Process Mixture of skew Normals}
\usage{
DPM_GibbsSampler_SkewN(z, hyperG0, a, b, N, doPlot = TRUE,
  nbclust_init = 30, plotevery = 1, ...)
}
\arguments{
  \item{z}{data matrix \code{d x n} with \code{d}
  dimensions in rows and \code{n} observations in columns.}

  \item{hyperG0}{prior mixing distribution.}

  \item{a}{shape hyperparameter of the Gamma prior on the
  parameter of the Dirichlet Process.}

  \item{b}{scale hyperparameter of the Gamma prior on the
  parameter of the Dirichlet Process.}

  \item{N}{number of MCMC iterations.}

  \item{doPlot}{logical flag indicating wether to plot MCMC
  iteration or not. Default to \code{TRUE}.}

  \item{nbclust_init}{number of clusters at initialisation.
  Default to 30 (or less if there are less than 30
  observations).}
}
\description{
Hyperprior on the concentration parameter of the DP
}
\examples{
rm(list=ls())
library(ggplot2)
#Number of data
n <- 1000
#n <- 2000
set.seed(1234)
#set.seed(123)
#set.seed(4321)

# Sample data
xi <- matrix(nrow=2, ncol=4, c(-1.5, 1, 1.5, 1, 1.5, -2, -2, -2))
psi <- matrix(nrow=2, ncol=4, c(0.4, -0.6, 0.8, 0, 0.3, -0.7, -0.3, -1.2))
p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters

sdev <- array(dim=c(2,2,4))
sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
sdev[, ,4] <- .3*diag(2)


c <- rep(0,n)
z <- matrix(0, nrow=2, ncol=n)
for(k in 1:n){
 c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
 z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + sdev[, , c[k]]\%*\%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
 #cat(k, "/", n, " observations simulated\\n", sep="")
}

 d<-2
 # Set parameters of G0
 hyperG0 <- list()
 hyperG0[["b_xi"]] <- rep(0,d)
 hyperG0[["b_psi"]] <- rep(0,d)
 hyperG0[["kappa"]] <- 0.001
 hyperG0[["D_xi"]] <- 0.1
 hyperG0[["D_psi"]] <- 0.1
 hyperG0[["nu"]] <- d+2
 hyperG0[["lambda"]] <- diag(d)/10

 # hyperprior on the Scale parameter of DPM
 a <- 0.0001
 b <- 0.0001

 # do some plots
 doPlot <- TRUE
 nbclust_init <- 30



 ## Data
 ########
 p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
       + geom_point()
       + ggtitle("Simple example in 2d data")
       +xlab("D1")
       +ylab("D2")
       +theme_bw())
 p


 ## alpha priors plots
 #####################
 prioralpha <- data.frame("alpha"=rgamma(n=5000, shape=a, scale=1/b),
                         "distribution" =factor(rep("prior",5000),
                         levels=c("prior", "posterior")))
 p <- (ggplot(prioralpha, aes(x=alpha))
       + geom_histogram(aes(y=..density..),
                        colour="black", fill="white")
       + geom_density(alpha=.2, fill="red")
       + ggtitle(paste("Prior distribution on alpha: Gamma(", a,
                 ",", b, ")\\n", sep=""))
      )
 p



 # Gibbs sampler for Dirichlet Process Mixtures
 ##############################################
 MCMCsample_sn <- DPM_GibbsSampler_SkewN(z, hyperG0, a, b, N=30, doPlot, nbclust_init, plotevery=1, gg.add=list(theme_bw()))
 plot_ConvDPM(MCMCsample_sn, from=2)

 hyperG0[["mu"]] <- rep(0,d)
 MCMCsample_n <- gibbsDPMsliceprior(z, hyperG0, a, b, N=100, doPlot, nbclust_init, plotevery=5)
 plot_ConvDPM(MCMCsample_n, from=2)





 # k-means

 plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster,
      xlab = "d = 1", ylab= "d = 2", main="k-means with K=4 clusters")

 KM <- kmeans(t(z), centers=4)
 dataKM <- data.frame("X"=z[1,], "Y"=z[2,],
                    "Cluster"=as.character(KM$cluster))
 dataCenters <- data.frame("X"=KM$centers[,1],
                           "Y"=KM$centers[,2],
                           "Cluster"=rownames(KM$centers))

 p <- (ggplot(dataKM)
       + geom_point(aes(x=X, y=Y, col=Cluster))
       + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster),
                    data=dataCenters, shape=22, size=5)
       + scale_colour_discrete(name="Cluster")
       + ggtitle("K-means with K=4 clusters\\n"))
 p

 postalpha <- data.frame("alpha"=MCMCsample$alpha[50:500],
                         "distribution" = factor(rep("posterior",500-49),
                         levels=c("prior", "posterior")))
 p <- (ggplot(postalpha, aes(x=alpha))
       + geom_histogram(aes(y=..density..), binwidth=.1,
                        colour="black", fill="white")
       + geom_density(alpha=.2, fill="blue")
       + ggtitle("Posterior distribution of alpha\\n")
       + geom_vline(aes(xintercept=mean(alpha, na.rm=T)),   # Ignore NA values for mean
                    color="red", linetype="dashed", size=1)  # Overlay with transparent density plot
     )
 p

 p <- (ggplot(drop=FALSE, alpha=.6)
       + geom_density(aes(x=alpha, fill=distribution),
                      color=NA, alpha=.6,
                      data=prioralpha)
       + geom_density(aes(x=alpha, fill=distribution),
                      color=NA, alpha=.6,
                      data=postalpha)
       + ggtitle("Prior and posterior distributions of alpha\\n")
       + scale_fill_discrete(drop=FALSE)
     )
 p
}
\author{
Boris Hejblum
}

