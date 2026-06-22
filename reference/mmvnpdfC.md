# C++ implementation of multivariate Normal probability density function for multiple inputs

C++ implementation of multivariate Normal probability density function
for multiple inputs

## Usage

``` r
mmvnpdfC(x, mean, varcovM, Log = TRUE)
```

## Arguments

- x:

  data matrix of dimension `p x n`, `p` being the dimension of the data
  and n the number of data points.

- mean:

  mean vectors matrix of dimension `p x K`, `K` being the number of
  distributions for which the density probability has to be evaluated.

- varcovM:

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension `K x n`.

## Examples

``` r
if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
               mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
               mmvnpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), Log=FALSE),
               times=1000L)
microbenchmark(mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
                      varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
               mvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), mean=c(-0.2, 0.3),
                       varcovM=matrix(c(2, 0.2, 0.2, 2), ncol=2), Log=FALSE),
               mmvnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
                        mean=matrix(c(-0.2, 0.3), nrow=2, ncol=1),
                        varcovM=list(matrix(c(2, 0.2, 0.2, 2), ncol=2)), Log=FALSE),
               times=1000L)
microbenchmark(mvnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
                      mean=list(c(0,0),c(-1,-1), c(1.5,1.5)),
                      varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
               mmvnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
                        mean=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
                        varcovM=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
               times=1000L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Loading required package: microbenchmark
#> Unit: microseconds
#>                                                                                                                                                                                                        expr
#>            mvnpdf(x = matrix(c(rep(1.96, 2), rep(0, 2)), nrow = 2, ncol = 2),      mean = list(c(0, 0), c(-1, -1), c(1.5, 1.5)), varcovM = list(diag(2),          10 * diag(2), 20 * diag(2)), Log = FALSE)
#>  mmvnpdfC(matrix(c(rep(1.96, 2), rep(0, 2)), nrow = 2, ncol = 2),      mean = matrix(c(0, 0, -1, -1, 1.5, 1.5), nrow = 2, ncol = 3),      varcovM = list(diag(2), 10 * diag(2), 20 * diag(2)), Log = FALSE)
#>      min       lq      mean  median       uq      max neval
#>  218.067 225.9165 241.96054 232.494 245.9845 3289.780  1000
#>   11.391  12.5235  15.72114  16.411  17.7430   67.075  1000
```
