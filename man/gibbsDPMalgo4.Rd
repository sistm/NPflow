\name{gibbsDPMalgo4}
\alias{gibbsDPMalgo4}
\title{Gibbs Sapling DPM Alghorithm 4: slice sampling}
\usage{
gibbsDPMalgo4(z, hyperG0, alpha, N, doPlot = TRUE)
}
\arguments{
  \item{z}{data}

  \item{hyperG0}{prior mixing distribution}

  \item{alpha}{parameter of the Dirichlet Process}

  \item{N}{number of MCMC iterations}
}
\description{
Gibbs Sapling DPM Alghorithm 4: slice sampling
}
\examples{
rm(list=ls())
#Number of data
n <- 100
#set.seed(1231)

# Sample data
m <- matrix(nrow=2, ncol=4, c(-1, 1, 0, 2, 1, -2, -1, -2))
p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters

library(expm)
s <- array(dim=c(2,2,4))
s[, ,1] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.1))
s[, ,2] <- matrix(nrow=2, ncol=2, c(0.01, 0, 0, 0.1))
s[, ,3] <- matrix(nrow=2, ncol=2, c(0.1, 0.08, 0.08, 0.1))
s[, ,4] <- .1*diag(2)
c <- rep(0,n)
z <- matrix(0, nrow=2, ncol=n)
for(k in 1:n){
 c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
 z[,k] <- m[, c[k]] + sqrtm(s[, , c[k]])\%*\%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
}

 # Set parameters of G0
 hyperG0 <- list()
 hyperG0[["mu"]] <- c(0,0)
 hyperG0[["kappa"]] <- 1
 hyperG0[["nu"]] <- 4
 hyperG0[["lambda"]] <- diag(2)
 # Scale parameter of DPM
 alpha <- 4
 # Number of iterations
 N <- 10
 # do some plots
 doPlot <- TRUE

 # Gibbs sampler for Dirichlet Process Mixtures
 GibSample <- gibbsDPMalgo4(z, hyperG0, alpha, N, doPlot)
}
\author{
Boris Hejblum
}

