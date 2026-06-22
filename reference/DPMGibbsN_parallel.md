# Slice Sampling of the Dirichlet Process Mixture Model with a prior on alpha

Slice Sampling of the Dirichlet Process Mixture Model with a prior on
alpha

## Usage

``` r
DPMGibbsN_parallel(
  Ncpus,
  type_connec,
  z,
  hyperG0,
  a = 1e-04,
  b = 1e-04,
  N,
  doPlot = TRUE,
  nbclust_init = 30,
  plotevery = N/10,
  diagVar = TRUE,
  use_variance_hyperprior = TRUE,
  verbose = TRUE,
  monitorfile = "",
  ...
)
```

## Arguments

- Ncpus:

  the number of processors available

- type_connec:

  The type of connection between the processors. Supported cluster types
  are `"SOCK"`, `"FORK"`, `"MPI"`, and `"NWS"`. See also
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html).

- z:

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

- use_variance_hyperprior:

  logical flag indicating whether a hyperprior is added for the variance
  parameter. Default is `TRUE` which decrease the impact of the variance
  prior on the posterior. `FALSE` is useful for using an informative
  prior.

- verbose:

  logical flag indicating whether partition info is written in the
  console at each MCMC iteration.

- monitorfile:

  a writable [connections](https://rdrr.io/r/base/connections.html) or a
  character string naming a file to write into, to monitor the progress
  of the analysis. Default is `""` which is no monitoring. See Details.

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

- `listU_mu`::

  a list of length `N` containing the matrices of mean vectors for all
  the mixture components at each MCMC iteration

- `listU_Sigma`::

  a list of length `N` containing the arrays of covariances matrices for
  all the mixture components at each MCMC iteration

- `U_SS_list`::

  a list of length `N` containing the lists of sufficient statistics for
  all the mixture components at each MCMC iteration

- `weights_list`::

  a list of length `N` containing the logposterior values at each MCMC
  iterations

- `logposterior_list`::

  a list of length `N` containing the logposterior values at each MCMC
  iterations

- `data`::

  the data matrix `d x n` with `d` dimensions in rows and `n`
  observations in columns

- `nb_mcmcit`::

  the number of MCMC iterations

- `clust_distrib`::

  the parametric distribution of the mixture component - `"gaussian"`

- `hyperG0`::

  the prior on the cluster location

## See also

[`DPMGibbsN`](http://sistm.github.io/NPflow/reference/DPMGibbsN.md)

## Author

Boris Hejblum

## Examples

``` r

# Scaling up: ----
rm(list=ls())
#Number of data
n <- 2000
set.seed(1234)

# Sample data
d <- 3
nclust <- 5
m <- matrix(nrow=d, ncol=nclust, runif(d*nclust)*8)
# p: cluster probabilities
p <- runif(nclust)
p <- p/sum(p)

# Covariance matrix of the clusters
sdev <- array(dim=c(d, d, nclust))
for (j in 1:nclust){
    sdev[, ,j] <- matrix(NA, nrow=d, ncol=d)
    diag(sdev[, ,j]) <- abs(rnorm(n=d, mean=0.3, sd=0.1))
    sdev[, ,j][lower.tri(sdev[, ,j], diag = FALSE)] <- rnorm(n=d*(d-1)/2,
    mean=0, sd=0.05)
    sdev[, ,j][upper.tri(sdev[, ,j], diag = FALSE)] <- (sdev[, ,j][
                                                        lower.tri(sdev[, ,j], diag = FALSE)])
}
c <- rep(0,n)
z <- matrix(0, nrow=d, ncol=n)
for(k in 1:n){
    c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
    z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
    #cat(k, "/", n, " observations simulated\n", sep="")
}

# hyperprior on the Scale parameter of DPM
a <- 0.001
b <- 0.001

# Number of iterations
N <- 25

# do some plots
doPlot <- TRUE

# Set parameters of G0
hyperG0 <- list()
hyperG0[["mu"]] <- rep(0, d)
hyperG0[["kappa"]] <- 0.01
hyperG0[["nu"]] <- d + 2
hyperG0[["lambda"]] <- diag(d)/10


nbclust_init <- 30

if(interactive()){
 library(doParallel)
 MCMCsample <- DPMGibbsN_parallel(Ncpus=2, type_connec="FORK", z, hyperG0, a, b, 
                                  N=1000, doPlot=FALSE, nbclust_init=30, 
                                  plotevery=100, gg.add=list(ggplot2::theme_bw(),
                                  ggplot2::guides(shape = 
                                    ggplot2::guide_legend(override.aes = list(fill="grey45")))),
                                  diagVar=FALSE)
}

```
