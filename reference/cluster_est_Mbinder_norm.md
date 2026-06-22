# Point estimate of the partition using a modified Binder loss function

Get a point estimate of the partition using a modified Binder loss
function for Gaussian components

## Usage

``` r
cluster_est_Mbinder_norm(c, Mu, Sigma, lambda = 0, a = 1, b = a, logposterior)
```

## Arguments

- c:

  a list of vector of length `n`. `c[[j]][i]` is the cluster allocation
  of observation `i=1...n` at iteration `j=1...N`.

- Mu:

  is a list of length `n` composed of `p x l` matrices. Where `l` is the
  maximum number of components per partition.

- Sigma:

  is list of length `n` composed of arrays containing a maximum of `l`
  `p x p` covariance matrices.

- lambda:

  is a nonnegative tunning parameter allowing further control over the
  distance function. Default is 0.

- a:

  nonnegative constant seen as the unit cost for pairwise
  misclassification. Default is 1.

- b:

  nonnegative constant seen as the unit cost for the other kind of
  pairwise misclassification. Default is 1.

- logposterior:

  vector of logposterior corresponding to each partition from `c` used
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

## Details

Note that he current implementation only allows Gaussian components.

The modified Binder loss function takes into account the distance
between mixture components using \#'the Bhattacharyya distance.

## References

JW Lau, PJ Green, Bayesian Model-Based Clustering Procedures, *Journal
of Computational and Graphical Statistics*, 16(3):526-558, 2007.

DA Binder, Bayesian cluster analysis, *Biometrika* 65(1):31-38, 1978.

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)
[`similarityMatC`](http://sistm.github.io/NPflow/reference/similarityMatC.md)
[`similarityMat_nocostC`](http://sistm.github.io/NPflow/reference/similarityMat_nocostC.md)

## Author

Chariff Alkhassim
