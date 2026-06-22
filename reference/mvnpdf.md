# multivariate-Normal probability density function

multivariate-Normal probability density function

## Usage

``` r
mvnpdf(x, mean, varcovM, Log = TRUE)
```

## Arguments

- x:

  p x n data matrix with n the number of observations and p the number
  of dimensions

- mean:

  mean vector or list of mean vectors (either a vector, a matrix or a
  list)

- varcovM:

  variance-covariance matrix or list of variance-covariance matrices
  (either a matrix or a list)

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## See also

`mvnpdf`,
[`mmvnpdfC`](http://sistm.github.io/NPflow/reference/mmvnpdfC.md)

## Author

Boris P. Hejblum

## Examples

``` r

mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)
#> [1] 0.05844094
dnorm(1.96)
#> [1] 0.05844094

mvnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      mean=c(0, 0), varcovM=diag(2), Log=FALSE
)
#> [1] 0.003415344
```
