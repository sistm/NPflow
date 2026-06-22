# Posterior estimation for Dirichlet process mixture of multivariate (potentially skew) distributions models

Partially collapse slice Gibbs sampling for Dirichlet process mixture of
multivariate normal, skew normal or skew t distributions.

## Usage

``` r
DPMpost(
  data,
  hyperG0,
  a = 1e-04,
  b = 1e-04,
  N,
  doPlot = TRUE,
  nbclust_init = 30,
  plotevery = floor(N/10),
  diagVar = TRUE,
  verbose = TRUE,
  distrib = c("gaussian", "skewnorm", "skewt"),
  ncores = 1,
  type_connec = "SOCK",
  informPrior = NULL,
  ...
)
```

## Arguments

- data:

  data matrix `d x n` with `d` dimensions in rows and `n` observations
  in columns.

- hyperG0:

  prior mixing distribution.

- a:

  shape hyperparameter of the Gamma prior on the concentration parameter
  of the Dirichlet Process. Default is `0.0001`.

- b:

  scale hyperparameter of the Gamma prior on the concentration parameter
  of the Dirichlet Process. Default is `0.0001`. If `0`, then the
  concentration is fixed set to `a`.

- N:

  number of MCMC iterations.

- doPlot:

  logical flag indicating whether to plot MCMC iteration or not. Default
  to `TRUE`.

- nbclust_init:

  number of clusters at initialization. Default to 30 (or less if there
  are less than 30 observations).

- plotevery:

  an integer indicating the interval between plotted iterations when
  `doPlot` is `TRUE`.

- diagVar:

  logical flag indicating whether the variance of each cluster is
  estimated as a diagonal matrix, or as a full matrix. Default is `TRUE`
  (diagonal variance).

- verbose:

  logical flag indicating whether partition info is written in the
  console at each MCMC iteration.

- distrib:

  the distribution used for the clustering. Current possibilities are
  `"gaussian"`, `"skewnorm"` and `"skewt"`.

- ncores:

  number of cores to use.

- type_connec:

  The type of connection between the processors. Supported cluster types
  are `"SOCK"`, `"FORK"`, `"MPI"`, and `"NWS"`. See also
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html).

- informPrior:

  an optional informative prior such as the approximation computed by
  `summary.DPMMclust`.

- ...:

  additional arguments to be passed to
  [`plot_DPM`](http://sistm.github.io/NPflow/reference/plot_DPM.md).
  Only used if `doPlot` is `TRUE`.

## Value

a object of class `DPMclust` with the following attributes:

- `mcmc_partitions`::

  a list of length `N`. Each element `mcmc_partitions[n]` is a vector of
  length `n` giving the partition of the `n` observations.

- `alpha`::

  a vector of length `N`. `cost[j]` is the cost associated to partition
  `c[[j]]`

- `U_SS_list`::

  a list of length `N` containing the lists of sufficient statistics for
  all the mixture components at each MCMC iteration

- `weights_list`::

  a list of length `N` containing the weights of each mixture component
  for each MCMC iterations

- `logposterior_list`::

  a list of length `N` containing the logposterior values at each MCMC
  iterations

- `data`::

  the data matrix `d x n` with `d` dimensions in rows and `n`
  observations in columns

- `nb_mcmcit`::

  the number of MCMC iterations

- `clust_distrib`::

  the parametric distribution of the mixture component

- `hyperG0`::

  the prior on the cluster location

## Details

This function is a wrapper around the following functions:
[`DPMGibbsN`](http://sistm.github.io/NPflow/reference/DPMGibbsN.md),
[`DPMGibbsN_parallel`](http://sistm.github.io/NPflow/reference/DPMGibbsN_parallel.md),
[`DPMGibbsN_SeqPrior`](http://sistm.github.io/NPflow/reference/DPMGibbsN_SeqPrior.md),
[`DPMGibbsSkewN`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewN.md),
[`DPMGibbsSkewN_parallel`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewN_parallel.md),
[`DPMGibbsSkewT`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewT.md),
[`DPMGibbsSkewT_parallel`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewT_parallel.md),
[`DPMGibbsSkewT_SeqPrior`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewT_SeqPrior.md),
[`DPMGibbsSkewT_SeqPrior_parallel`](http://sistm.github.io/NPflow/reference/DPMGibbsSkewT_SeqPrior_parallel.md).

## References

Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019)
Sequential Dirichlet Process Mixtures of Multivariate Skew
t-distributions for Model-based Clustering of Flow Cytometry Data. The
Annals of Applied Statistics, 13(1): 638-660. \<doi:
10.1214/18-AOAS1209\> \<arXiv: 1702.04407\>
<https://arxiv.org/abs/1702.04407>
[doi:10.1214/18-AOAS1209](https://doi.org/10.1214/18-AOAS1209)

## See also

[`summary.DPMMclust`](http://sistm.github.io/NPflow/reference/summary.DPMMclust.md)

## Author

Boris Hejblum

## Examples

``` r
#rm(list=ls())
set.seed(123)

# Exemple in 2 dimensions with skew-t distributions

# Generate data:
n <- 2000 # number of data points
d <- 2 # dimensions
ncl <- 4 # number of true clusters
sdev <- array(dim=c(d,d,ncl))
xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1.5, 1.5, 1.5, 2, -2.5, -2.5, -3))
psi <- matrix(nrow=d, ncol=4, c(0.3, -0.7, -0.8, 0, 0.3, -0.7, 0.2, 0.9))
nu <- c(100,25,8,5)
proba <- c(0.15, 0.05, 0.5, 0.3) # cluster frequencies
sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.2))
sdev[, ,4] <- .3*diag(2)
c <- rep(0,n)
w <- rep(1,n)
z <- matrix(0, nrow=d, ncol=n)
for(k in 1:n){
 c[k] = which(rmultinom(n=1, size=1, prob=proba)!=0)
 w[k] <- rgamma(1, shape=nu[c[k]]/2, rate=nu[c[k]]/2)
 z[,k] <- xi[, c[k]] + psi[, c[k]]*rtruncnorm(n=1, a=0, b=Inf, mean=0, sd=1/sqrt(w[k])) +
               (sdev[, , c[k]]/sqrt(w[k]))%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
}

# Define hyperprior
hyperG0 <- list()
hyperG0[["b_xi"]] <- rowMeans(z)
hyperG0[["b_psi"]] <- rep(0,d)
hyperG0[["kappa"]] <- 0.001
hyperG0[["D_xi"]] <- 100
hyperG0[["D_psi"]] <- 100
hyperG0[["nu"]] <- d+1
hyperG0[["lambda"]] <- diag(apply(z,MARGIN=1, FUN=var))/3


if(interactive()){
# Plot data
cytoScatter(z)

# Estimate posterior
MCMCsample_st <- DPMpost(data=z, hyperG0=hyperG0, N=2000,
   distrib="skewt",
   gg.add=list(ggplot2::theme_bw(),
      ggplot2::guides(shape=ggplot2::guide_legend(override.aes = list(fill="grey45"))))
 )
 s <- summary(MCMCsample_st, burnin = 1600, thin=5, lossFn = "Binder")
 s
 plot(s)
 #plot(s, hm=TRUE) # this can take a few sec...
 
 
 # more data plotting:
 library(ggplot2)
 p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
       + geom_point()
       + ggtitle("Unsupervised data")
       + xlab("D1")
       + ylab("D2")
       + theme_bw()
 )
 p

 c2plot <- factor(c)
 levels(c2plot) <- c("4", "1", "3", "2")
 pp <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,], "Cluster"=as.character(c2plot)))
       + geom_point(aes(x=X, y=Y, colour=Cluster, fill=Cluster))
       + ggtitle("True clusters")
       + xlab("D1")
       + ylab("D2")
       + theme_bw()
       + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6, shape=22)))
 )
 pp
}




# Exemple in 2 dimensions with Gaussian distributions

set.seed(1234)

# Generate data 
n <- 2000 # number of data points
d <- 2 # dimensions
ncl <- 4 # number of true clusters
m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2)) # cluster means
sdev <- array(dim=c(2, 2, 4)) # cluster standard-deviations
sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
sdev[, ,4] <- .3*diag(2)
proba <- c(0.15, 0.05, 0.5, 0.3) # cluster frequencies
c <- rep(0,n)
z <- matrix(0, nrow=2, ncol=n)
for(k in 1:n){
 c[k] = which(rmultinom(n=1, size=1, prob=proba)!=0)
 z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
}

# Define hyperprior
hyperG0 <- list()
hyperG0[["mu"]] <- rep(0,d)
hyperG0[["kappa"]] <- 0.001
hyperG0[["nu"]] <- d+2
hyperG0[["lambda"]] <- diag(d)


if(interactive()){
# Plot data
cytoScatter(z)

# Estimate posterior
MCMCsample_n <- DPMpost(data=z, hyperG0=hyperG0, N=2000,
   distrib="gaussian", diagVar=FALSE,
   gg.add=list(ggplot2::theme_bw(),
      ggplot2::guides(shape=ggplot2::guide_legend(override.aes = list(fill="grey45"))))
 )
 s <- summary(MCMCsample_n, burnin = 1500, thin=5, lossFn = "Binder")
 s
 plot(s)
 #plot(s, hm=TRUE) # this can take a few sec...
 
 
 # more data plotting:
 library(ggplot2)
 p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
       + geom_point()
       + ggtitle("Unsupervised data")
       + xlab("D1")
       + ylab("D2")
       + theme_bw()
 )
 p

 c2plot <- factor(c)
 levels(c2plot) <- c("4", "1", "3", "2")
 pp <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,], "Cluster"=as.character(c2plot)))
       + geom_point(aes(x=X, y=Y, colour=Cluster, fill=Cluster))
       #+ ggtitle("Slightly overlapping skew-normal simulation\n")
       + xlab("D1")
       + ylab("D2")
       + theme_bw()
       + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6, shape=22)))
       + ggtitle("True clusters")
 )
 pp
}

```
