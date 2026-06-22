# Point estimate of the partition for the Binder loss function

Get a point estimate of the partition using the Binder loss function.

## Usage

``` r
cluster_est_binder(c, logposterior)
```

## Arguments

- c:

  a list of vector of length `n`. `c[[j]][i]` is the cluster allocation
  of observation `i=1...n` at iteration `j=1...N`.

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

## References

F Caron, YW Teh, TB Murphy, Bayesian nonparametric Plackett-Luce models
for the analysis of preferences for college degree programmes, *Annals
of Applied Statistics*, 8(2):1145-1181, 2014.

DB Dahl, Model-Based Clustering for Expression Data via a Dirichlet
Process Mixture Model, *Bayesian Inference for Gene Expression and
Proteomics*, K-A Do, P Muller, M Vannucci (Eds.), Cambridge University
Press, 2006.

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)
[`similarityMatC`](http://sistm.github.io/NPflow/reference/similarityMatC.md)

## Author

Francois Caron, Boris Hejblum
