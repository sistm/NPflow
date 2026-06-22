# C++ implementation of multivariate Normal probability density function for multiple inputs

C++ implementation of multivariate Normal probability density function
for multiple inputs

## Usage

``` r
mvnlikC(x, c, clustval, mu, sigma, loglik = TRUE)
```

## Arguments

- x:

  data matrix of dimension p x n, p being the dimension of the data and
  n the number of data points

- c:

  integer vector of cluster allocations with values from 1 to K

- clustval:

  vector of unique values from c in the order corresponding to the
  storage of cluster parameters in `xi`, `psi`, and `varcovM`

- mu:

  mean vectors matrix of dimension p x K, K being the number of clusters

- sigma:

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- loglik:

  logical flag or returning the log-likelihood instead of the
  likelihood. Default is `TRUE`.

## Value

a list:

- `"indiv"`::

  vector of likelihood of length n;

- `"clust"`::

  vector of likelihood of length K;

- `"total"`::

  total (log)-likelihood;

## Author

Boris Hejblum
