# Computes the co-clustering (or similarity) matrix

Computes the co-clustering (or similarity) matrix

## Usage

``` r
similarityMat(c, step = 1)
```

## Arguments

- c:

  a list of vector of length `n`. `c[[j]][i]` is the cluster allocation
  of observation `i=1...n` at iteration `j=1...N`.

- step:

  provide co-clustering every `step` iterations. Default is 1.

## Value

A matrix of size `n x n` whose term `[i,j]` is the proportion of MCMC
iterations where observation `i` and observations `j` are allocated to
the same cluster.

## Author

Boris Hejblum
