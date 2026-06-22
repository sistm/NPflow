# Burning MCMC iterations from a Dirichlet Process Mixture Model.

Utility function for burning MCMC iteration from a `DPMMclust` object.

## Usage

``` r
burn.DPMMclust(x, burnin = 0, thin = 1)
```

## Arguments

- x:

  a `DPMMclust` object.

- burnin:

  the number of MCMC iterations to burn (default is `0`).

- thin:

  the spacing at which MCMC iterations are kept. Default is `1`, i.e. no
  thining.

## Value

a `DPMMclust` object minus the burnt iterations

## See also

[`summary.DPMMclust`](http://sistm.github.io/NPflow/reference/summary.DPMMclust.md)

## Author

Boris Hejblum
