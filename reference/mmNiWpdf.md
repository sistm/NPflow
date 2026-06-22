# multivariate Normal inverse Wishart probability density function for multiple inputs

multivariate Normal inverse Wishart probability density function for
multiple inputs

## Usage

``` r
mmNiWpdf(mu, Sigma, U_mu0, U_kappa0, U_nu0, U_lambda0, Log = TRUE)
```

## Arguments

- mu:

  data matrix of dimension `p x n`, `p` being the dimension of the data
  and n the number of data points, where each column is an observed mean
  vector.

- Sigma:

  list of length `n` of observed variance-covariance matrices, each of
  dimensions `p x p`.

- U_mu0:

  mean vectors matrix of dimension `p x K`, `K` being the number of
  distributions for which the density probability has to be evaluated

- U_kappa0:

  vector of length `K` of scale parameters.

- U_nu0:

  vector of length `K` of degree of freedom parameters.

- U_lambda0:

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension K x n
