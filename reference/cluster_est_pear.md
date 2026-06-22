# Gets a point estimate of the partition using posterior expected adjusted Rand index (PEAR)

Gets a point estimate of the partition using posterior expected adjusted
Rand index (PEAR)

## Usage

``` r
cluster_est_pear(c)
```

## Arguments

- c:

  a list of vector of length `n`. `c[[j]][i]` is the cluster allocation
  of observation `i=1...n` at iteration `j=1...N`.

## Value

a `list`:

- `c_est`::

  a vector of length `n`. Point estimate of the partition

- `pear`::

  a vector of length `N`. `pear[j]` is the posterior expected adjusted
  Rand index associated to partition `c[[j]]`

- `similarity`::

  matrix of size `n x n`. Similarity matrix (see
  [`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md))

- `opt_ind`::

  the index of the optimal partition among the MCMC iterations.

## References

A. Fritsch, K. Ickstadt. Improved Criteria for Clustering Based on the
Posterior Similarity Matrix, in Bayesian Analysis, Vol.4 : p.367-392
(2009)

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)
[`similarityMatC`](http://sistm.github.io/NPflow/reference/similarityMatC.md)

## Author

Chariff Alkhassim
