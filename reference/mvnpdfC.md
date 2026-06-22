# C++ implementation of multivariate normal probability density function for multiple inputs

Based on the implementation from Nino Hardt and Dicko Ahmadou
<https://gallery.rcpp.org/articles/dmvnorm_arma/> (accessed in August
2014)

## Usage

``` r
mvnpdfC(x, mean, varcovM, Log = TRUE)
```

## Arguments

- x:

  data matrix

- mean:

  mean vector

- varcovM:

  variance covariance matrix

- Log:

  logical flag for returning the log of the probability density
  function. Default is `TRUE`

## Value

vector of densities

## Author

Boris P. Hejblum

## Examples

``` r
mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#> [1] 0.05844094
mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#> [1] 0.05844094
mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#> [1] -2.839739
mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1))
#> [1] -2.839739

if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(dnorm(1.96),
               mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
               mvnpdfC(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE),
               times=10000L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Unit: nanoseconds
#>                                                                 expr   min
#>                                                          dnorm(1.96)   732
#>   mvnpdf(x = matrix(1.96), mean = 0, varcovM = diag(1), Log = FALSE) 48962
#>  mvnpdfC(x = matrix(1.96), mean = 0, varcovM = diag(1), Log = FALSE)  4198
#>     lq      mean median    uq     max neval
#>    922  1296.116   1302  1542  109715 10000
#>  54552 58665.492  55644 57507 3175136 10000
#>   4839  6225.256   5961  6523 2888742 10000
```
