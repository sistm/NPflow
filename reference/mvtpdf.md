# multivariate Student's t-distribution probability density function

multivariate Student's t-distribution probability density function

## Usage

``` r
mvtpdf(x, mean, varcovM, df, Log = TRUE)
```

## Arguments

- x:

  `p x n` data matrix with `n` the number of observations and `p` the
  number of dimensions

- mean:

  mean vector or list of mean vectors (either a vector, a matrix or a
  list)

- varcovM:

  variance-covariance matrix or list of variance-covariance matrices
  (either a matrix or a list)

- df:

  a numeric vector or a list of the degrees of freedom (either a vector
  or a list)

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Examples

``` r
mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000)
#> [1] -2.839738
mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#> [1] -2.839739

mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10)
#> [1] -2.731911

mvtpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      mean=c(0, 0), varcovM=diag(2), df=10
)
#> [1] -5.258057
```
