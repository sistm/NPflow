# C++ implementation of multivariate skew Normal probability density function for multiple inputs

C++ implementation of multivariate skew Normal probability density
function for multiple inputs

## Usage

``` r
mmvsnpdfC(x, xi, psi, sigma, Log = TRUE)
```

## Arguments

- x:

  data matrix of dimension `p x n`, `p` being the dimension of the data
  and n the number of data points.

- xi:

  mean vectors matrix of dimension `p x K`, `K` being the number of
  distributions for which the density probability has to be evaluated.

- psi:

  skew parameter vectors matrix of dimension `p x K`.

- sigma:

  list of length K of variance-covariance matrices, each of dimensions
  `p x p`.

- Log:

  logical flag for returning the log of the probability density
  function. Default is `TRUE`.

## Value

matrix of densities of dimension `K x n`.

## Author

Boris Hejblum

## Examples

``` r
mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)), Log=FALSE
         )
#>            [,1]
#> [1,] 0.05046623
mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2))
         )
#>           [,1]
#> [1,] -2.986451

if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), xi=c(0, 0), psi=c(1, 1),
                       sigma=diag(2), Log=FALSE),
               mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1), xi=matrix(c(0, 0)),
                         psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)), Log=FALSE),
               times=1000L
             )
microbenchmark(mvsnpdf(x=matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
                      xi=list(c(0,0),c(-1,-1), c(1.5,1.5)),
                      psi=list(c(0.1,0.1),c(-0.1,-1), c(0.5,-1.5)),
                      sigma=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
               mmvsnpdfC(matrix(c(rep(1.96,2),rep(0,2)), nrow=2, ncol=2),
                         xi=matrix(c(0,0,-1,-1, 1.5,1.5), nrow=2, ncol=3),
                         psi=matrix(c(0.1,0.1,-0.1,-1, 0.5,-1.5), nrow=2, ncol=3),
                         sigma=list(diag(2),10*diag(2), 20*diag(2)), Log=FALSE),
              times=1000L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Unit: microseconds
#>                                                                                                                                                                                                                                                                                       expr
#>                                 mvsnpdf(x = matrix(c(rep(1.96, 2), rep(0, 2)), nrow = 2, ncol = 2),      xi = list(c(0, 0), c(-1, -1), c(1.5, 1.5)), psi = list(c(0.1,          0.1), c(-0.1, -1), c(0.5, -1.5)), sigma = list(diag(2),          10 * diag(2), 20 * diag(2)), Log = FALSE)
#>  mmvsnpdfC(matrix(c(rep(1.96, 2), rep(0, 2)), nrow = 2, ncol = 2),      xi = matrix(c(0, 0, -1, -1, 1.5, 1.5), nrow = 2, ncol = 3),      psi = matrix(c(0.1, 0.1, -0.1, -1, 0.5, -1.5), nrow = 2,          ncol = 3), sigma = list(diag(2), 10 * diag(2), 20 * diag(2)),      Log = FALSE)
#>      min     lq      mean   median       uq      max neval
#>  358.910 372.22 395.76169 384.5480 397.4865 3367.586  1000
#>   14.697  16.05  19.95065  20.9745  22.3520   61.154  1000
```
