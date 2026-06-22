# multivariate skew-t probability density function

multivariate skew-t probability density function

## Usage

``` r
mvstpdf(x, xi, sigma, psi, df, Log = TRUE)
```

## Arguments

- x:

  `p x n` data matrix with `n` the number of observations and `p` the
  number of dimensions

- xi:

  mean vector or list of mean vectors (either a vector, a matrix or a
  list)

- sigma:

  variance-covariance matrix or list of variance-covariance matrices
  (either a matrix or a list)

- psi:

  skew parameter vector or list of skew parameter vectors (either a
  vector, a matrix or a list)

- df:

  a numeric vector or a list of the degrees of freedom (either a vector
  or a list)

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## See also

[`mvtpdf`](http://sistm.github.io/NPflow/reference/mvtpdf.md),
[`mvsnpdf`](http://sistm.github.io/NPflow/reference/mvsnpdf.md),
[`mmvstpdfC`](http://sistm.github.io/NPflow/reference/mmvstpdfC.md),
[`mvstlikC`](http://sistm.github.io/NPflow/reference/mvstlikC.md)

## Examples

``` r
mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
      df=100000000, Log=FALSE
)
#>            [,1]
#> [1,] 0.05046623
mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
      Log=FALSE
)
#>            [,1]
#> [1,] 0.05046623
mvstpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      xi=c(0, 0), psi=c(1, 1), sigma=diag(2),
      df=100000000
)
#>           [,1]
#> [1,] -2.986451
mvsnpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1),
      xi=c(0, 0), psi=c(1, 1), sigma=diag(2)
)
#>           [,1]
#> [1,] -2.986451

```
