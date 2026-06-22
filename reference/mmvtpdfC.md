# C++ implementation of multivariate Normal probability density function for multiple inputs

C++ implementation of multivariate Normal probability density function
for multiple inputs

## Usage

``` r
mmvtpdfC(x, mean, varcovM, df, Log = TRUE)
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

- df:

  vector of length `K` of degree of freedom parameters.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension `K x n`.

## Author

Boris Hejblum

## Examples

``` r
mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#> [1] 0.05844094
mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000, Log=FALSE)
#> [1] 0.05844095
mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10000000, Log=FALSE)
#>            [,1]
#> [1,] 0.05844095

mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#> [1] -2.839739
mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000)
#> [1] -2.839738
mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10000000)
#>           [,1]
#> [1,] -2.839738

mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10)
#> [1] -2.731911
mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)), df=10)
#>           [,1]
#> [1,] -2.731911


if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=1, Log=FALSE),
               mmvtpdfC(x=matrix(1.96), mean=matrix(0), varcovM=list(diag(1)),
                        df=c(1), Log=FALSE),
               times=10000L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Unit: microseconds
#>                                                                                                expr
#>                     mvtpdf(x = matrix(1.96), mean = 0, varcovM = diag(1), df = 1,      Log = FALSE)
#>  mmvtpdfC(x = matrix(1.96), mean = matrix(0), varcovM = list(diag(1)),      df = c(1), Log = FALSE)
#>     min      lq      mean  median     uq      max neval
#>  54.222 57.3070 61.751788 58.5845 60.343 3036.638 10000
#>   6.352  7.0935  8.448238  8.5260  9.258   29.205 10000
```
