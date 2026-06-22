# C++ implementation of multivariate Normal inverse Wishart probability density function for multiple inputs

C++ implementation of multivariate Normal inverse Wishart probability
density function for multiple inputs

## Usage

``` r
mmNiWpdfC(Mu, Sigma, U_Mu0, U_Kappa0, U_Nu0, U_Sigma0, Log = TRUE)
```

## Arguments

- Mu:

  data matrix of dimension `p x n`, `p` being the dimension of the data
  and n the number of data points, where each column is an observed mean
  vector.

- Sigma:

  list of length `n` of observed variance-covariance matrices, each of
  dimensions `p x p`.

- U_Mu0:

  mean vectors matrix of dimension `p x K`, `K` being the number of
  distributions for which the density probability has to be evaluated

- U_Kappa0:

  vector of length `K` of scale parameters.

- U_Nu0:

  vector of length `K` of degree of freedom parameters.

- U_Sigma0:

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension K x n

## References

Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019)
Sequential Dirichlet Process Mixtures of Multivariate Skew
t-distributions for Model-based Clustering of Flow Cytometry Data. The
Annals of Applied Statistics, 13(1): 638-660. \<doi:
10.1214/18-AOAS1209\>. \<arXiv: 1702.04407\>.
<https://arxiv.org/abs/1702.04407>
[doi:10.1214/18-AOAS1209](https://doi.org/10.1214/18-AOAS1209)
