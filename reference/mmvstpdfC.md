# C++ implementation of multivariate Normal probability density function for multiple inputs

C++ implementation of multivariate Normal probability density function
for multiple inputs

## Usage

``` r
mmvstpdfC(x, xi, psi, sigma, df, Log = TRUE)
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

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- df:

  vector of length K of degree of freedom parameters.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension `K x n`.

## Author

Boris Hejblum

## Examples

``` r
mmvstpdfC(x = matrix(c(3.399890,-5.936962), ncol=1), xi=matrix(c(0.2528859,-2.4234067)),
psi=matrix(c(11.20536,-12.51052), ncol=1),
sigma=list(matrix(c(0.2134011, -0.0382573, -0.0382573, 0.2660086), ncol=2)),
df=c(7.784106)
)
#>           [,1]
#> [1,] -3.207832
mvstpdf(x = matrix(c(3.399890,-5.936962), ncol=1), xi=c(0.2528859,-2.4234067),
psi=c(11.20536,-12.51052),
sigma=matrix(c(0.2134011, -0.0382573, -0.0382573, 0.2660086), ncol=2),
df=c(7.784106)
)
#>           [,1]
#> [1,] -3.207832

#skew-normal limit
mmvsnpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2))
         )
#>           [,1]
#> [1,] -2.986451
mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
       xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
       df=100000000
       )
#>           [,1]
#> [1,] -2.986451
mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1), sigma=list(diag(2)),
         df=100000000
         )
#>           [,1]
#> [1,] -2.986451

#non-skewed limit
mmvtpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
        mean=matrix(c(0, 0)), varcovM=list(diag(2)),
        df=10
        )
#>           [,1]
#> [1,] -5.258057
mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
         xi=matrix(c(0, 0)), psi=matrix(c(0, 0),ncol=1), sigma=list(diag(2)),
         df=10
         )
#>           [,1]
#> [1,] -5.258057

if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
                       xi=c(0, 0), psi=c(1, 1),
                       sigma=diag(2), df=10),
               mmvstpdfC(x=matrix(rep(1.96,2), nrow=2, ncol=1),
                         xi=matrix(c(0, 0)), psi=matrix(c(1, 1),ncol=1),
                         sigma=list(diag(2)), df=10),
               times=1000L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Unit: microseconds
#>                                                                                                                                                      expr
#>                                         mvstpdf(x = matrix(rep(1.96, 2), nrow = 2, ncol = 1), xi = c(0,      0), psi = c(1, 1), sigma = diag(2), df = 10)
#>  mmvstpdfC(x = matrix(rep(1.96, 2), nrow = 2, ncol = 1), xi = matrix(c(0,      0)), psi = matrix(c(1, 1), ncol = 1), sigma = list(diag(2)),      df = 10)
#>      min      lq      mean   median       uq      max neval
#>  134.912 138.859 150.27427 141.1485 147.4650 3079.377  1000
#>    9.488  10.330  13.16736  14.1915  15.0635   62.376  1000
```
