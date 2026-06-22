# Post-processing Dirichlet Process Mixture Models results to get a mixture distribution of the posterior locations

Post-processing Dirichlet Process Mixture Models results to get a
mixture distribution of the posterior locations

## Usage

``` r
postProcess.DPMMclust(
  x,
  burnin = 0,
  thin = 1,
  gs = NULL,
  lossFn = "F-measure",
  K = 10,
  ...
)
```

## Arguments

- x:

  a `DPMMclust` object.

- burnin:

  integer giving the number of MCMC iterations to burn (defaults is
  half)

- thin:

  integer giving the spacing at which MCMC iterations are kept. Default
  is `1`, i.e. no thining.

- gs:

  optional vector of length `n` containing the gold standard partition
  of the `n` observations to compare to the point estimate.

- lossFn:

  character string specifying the loss function to be used. Either
  "F-measure" or "Binder" (see Details). Default is "F-measure".

- K:

  integer giving the number of mixture components. Default is `10`.

- ...:

  further arguments passed to or from other methods

## Value

a `list`:

- `burnin`::

  an integer passing along the `burnin` argument

- `thin`::

  an integer passing along the `thin` argument

- `lossFn`::

  a character string passing along the `lossFn` argument

- `point_estim`::

- `loss`::

- `index_estim`::

## Details

The cost of a point estimate partition is calculated using either a
pairwise coincidence loss function (Binder), or 1-Fmeasure (F-measure).

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)
[`summary.DPMMclust`](http://sistm.github.io/NPflow/reference/summary.DPMMclust.md)

## Author

Boris Hejblum
