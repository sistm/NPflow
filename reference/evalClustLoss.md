# ELoss of a partition point estimate compared to a gold standard

Evaluate the loss of a point estimate of the partition compared to a
gold standard according to a given loss function

## Usage

``` r
evalClustLoss(c, gs, lossFn = "F-measure", a = 1, b = 1)
```

## Arguments

- c:

  vector of length `n` containing the estimated partition of the `n`
  observations.

- gs:

  vector of length `n` containing the gold standard partition of the `n`
  observations.

- lossFn:

  character string specifying the loss function to be used. Either
  "F-measure" or "Binder" (see Details). Default is "F-measure".

- a:

  only relevant if `lossFn` is "Binder". Penalty for wrong co-clustering
  in `c` compared to `gs`. Defaults is 1.

- b:

  only relevant if `lossFn` is "Binder". Penalty for missed
  co-clustering in `c` compared to `gs`. Defaults is 1.

## Value

the cost of the point estimate `c` in regard of the gold standard `gs`
for a given loss function.

## Details

The cost of a point estimate partition is calculated using either a
pairwise coincidence loss function (Binder), or 1-Fmeasure (F-measure).

## References

J.W. Lau & P.J. Green. Bayesian Model-Based Clustering Procedures,
Journal of Computational and Graphical Statistics, 16(3): 526-558, 2007.

D. B. Dahl. Model-Based Clustering for Expression Data via a Dirichlet
Process Mixture Model, in Bayesian Inference for Gene Expression and
Proteomics, K.-A. Do, P. Muller, M. Vannucci (Eds.), Cambridge
University Press, 2006.

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md),
[`cluster_est_binder`](http://sistm.github.io/NPflow/reference/cluster_est_binder.md)

## Author

Boris Hejblum
