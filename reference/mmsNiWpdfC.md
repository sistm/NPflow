# C++ implementation of multivariate structured Normal inverse Wishart probability density function for multiple inputs

C++ implementation of multivariate structured Normal inverse Wishart
probability density function for multiple inputs

## Usage

``` r
mmsNiWpdfC(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0, Log = TRUE)
```

## Arguments

- xi:

  data matrix of dimensions `p x n` where columns contain the observed
  mean vectors.

- psi:

  data matrix of dimensions `p x n` where columns contain the observed
  skew parameter vectors.

- Sigma:

  list of length `n` of observed variance-covariance matrices, each of
  dimensions `p x p`.

- U_xi0:

  mean vectors matrix of dimension `p x K`, `K` being the number of
  distributions for which the density probability has to be evaluated.

- U_psi0:

  skew parameter vectors matrix of dimension `p x K`.

- U_B0:

  list of length `K` of structured scale matrices, each of dimensions
  `p x p`.

- U_Sigma0:

  list of length `K` of variance-covariance matrices, each of dimensions
  `p x p`.

- U_df0:

  vector of length `K` of degree of freedom parameters.

- Log:

  logical flag for returning the log of the probability density
  function. Defaults is `TRUE`.

## Value

matrix of densities of dimension `K x n`

## References

Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019)
Sequential Dirichlet Process Mixtures of Multivariate Skew
t-distributions for Model-based Clustering of Flow Cytometry Data. The
Annals of Applied Statistics, 13(1): 638-660. \<doi:
10.1214/18-AOAS1209\>. \<arXiv: 1702.04407\>.
<https://arxiv.org/abs/1702.04407>
[doi:10.1214/18-AOAS1209](https://doi.org/10.1214/18-AOAS1209)
