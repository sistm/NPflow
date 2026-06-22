# Point estimate of the partition using the F-measure as the cost function.

Get a point estimate of the partition using the F-measure as the cost
function.

## Usage

``` r
cluster_est_Fmeasure(c, logposterior)
```

## Arguments

- c:

  a list of vector of length `n`. `c[[j]][i]` is the cluster allocation
  of observation `i=1...n` at iteration `j=1...N`.

- logposterior:

  a vector of logposterior corresponding to each partition from `c` used
  to break ties when minimizing the cost function

## Value

a `list`:

- `c_est`::

  a vector of length `n`. Point estimate of the partition

- `cost`::

  a vector of length `N`. `cost[j]` is the cost associated to partition
  `c[[j]]`

- `similarity`::

  matrix of size `n x n`. Similarity matrix (see
  [`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md))

- `opt_ind`::

  the index of the optimal partition among the MCMC iterations.

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)

## Author

Francois Caron, Boris Hejblum
