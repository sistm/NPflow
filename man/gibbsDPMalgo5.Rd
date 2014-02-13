\name{gibbsDPMalgo5}
\alias{gibbsDPMalgo5}
\title{Gibbs Sapling DPM Alghorithm 4: slice sampling
with a prior on alpha}
\usage{
gibbsDPMalgo5(z, hyperG0, a, b, N, doPlot = TRUE, nbclust_init = 30)
}
\arguments{
  \item{z}{data}

  \item{hyperG0}{prior mixing distribution}

  \item{a}{shape parameter f the Gamma hyperprior}

  \item{b}{scale parameter f the Gamma hyperprior}

  \item{N}{number of MCMC iterations}
}
\description{
Gibbs Sapling DPM Alghorithm 4: slice sampling with a prior
on alpha
}
\examples{
rm(list=ls())
#Number of data
n <- 2000
set.seed(1234)
#set.seed(4321)

# Sample data
m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2))
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
 cat(k, "/", n, " observations simulated\\n", sep="")
}

 # Set parameters of G0
 hyperG0 <- list()
 hyperG0[["mu"]] <- c(0,0)
 hyperG0[["kappa"]] <- 1
 hyperG0[["nu"]] <- 4
 hyperG0[["lambda"]] <- diag(2)
 # hyperprior on the Scale parameter of DPM
 a <- 0.001
 b <- 0.001
 plot(density(rgamma(n=5000, a, 1/b)))
 # Number of iterations
 N <- 15

 # do some plots
 doPlot <- TRUE
 nbclust_init <- 30

 # Gibbs sampler for Dirichlet Process Mixtures
 MCMCsample <- gibbsDPMalgo5(z, hyperG0, a, b, N, doPlot, nbclust_init)

 alpha_m <- mean(MCMCsample$alpha[floor(length(MCMCsample$alpha)/2):length(MCMCsample$alpha)])
 alpha_m*log(n/alpha_m) # order of the number of clusters
}
\author{
Boris Hejblum
}

