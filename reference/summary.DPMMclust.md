# Summarizing Dirichlet Process Mixture Models

Summary methods for `DPMMclust` objects.

## Usage

``` r
# S3 method for class 'DPMMclust'
summary(
  object,
  burnin = 0,
  thin = 1,
  gs = NULL,
  lossFn = "Binder",
  posterior_approx = FALSE,
  ...
)
```

## Arguments

- object:

  a `DPMMclust` object.

- burnin:

  integer giving the number of MCMC iterations to burn (defaults is
  half)

- thin:

  integer giving the spacing at which MCMC iterations are kept. Default
  is `1`, i.e. no thining.

- gs:

  optional vector of length `n` containing the gold standard partition
  of the `n` observations to compare to the point estimate

- lossFn:

  character string specifying the loss function to be used. Either
  "F-measure" or "Binder" (see Details). Default is "Binder".

- posterior_approx:

  logical flag whether a parametric approximation of the posterior
  should be computed. Default is `FALSE`

- ...:

  further arguments passed to or from other methods

## Value

a `list` containing the following elements:

- `nb_mcmcit`::

  an integer giving the value of `m`, the number of retained sampled
  partitions, i.e. `(N - burnin)/thin`

- `burnin`::

  an integer passing along the `burnin` argument

- `thin`::

  an integer passing along the `thin` argument

- `lossFn`::

  a character string passing along the `lossFn` argument

- `clust_distrib`::

  a character string passing along the `clust_distrib` argument

- `point_estim`::

  a `list` containing:

  `c_est`:

  :   a vector of length `n`containing the point estimated clustering
      for each observations

  `cost`:

  :   a vector of length `m` containing the cost of each sampled
      partition

  `Fmeas`:

  :   if `lossFn` is `'F-measure'`, the `m x m` matrix of total
      F-measures for each pair of sampled partitions

  `opt_ind`:

  :   the index of the point estimate partition among the `m` sampled

- `loss`::

  the loss for the point estimate. `NA` if `lossFn` is not `'Binder'`

- `param_posterior`::

  a list containing the parametric approximation of the posterior,
  suitable to be plugged in as prior for a new MCMC algorithm run

- `mcmc_partitions`::

  a list containing the `m` sampled partitions

- `alpha`::

  a vector of length `m` with the values of the `alpha` DP parameter

- `index_estim`::

  the index of the point estimate partition among the `m` sampled

- `hyperG0`::

  a list passing along the prior, i.e. the `hyperG0` argument

- `logposterior_list`::

  a list of length `m` containing the logposterior and its
  decomposition, for each sampled partition

- `U_SS_list`::

  a list of length `m` containing the containing the lists of sufficient
  statistics for all the mixture components, for each sampled partition

- `data`::

  a `d x n` matrix containing the clustered data

## Details

The cost of a point estimate partition is calculated using either a
pairwise coincidence loss function (Binder), or 1-Fmeasure (F-measure).

The number of retained sampled partitions is `m = (N - burnin)/thin`

## See also

[`similarityMat`](http://sistm.github.io/NPflow/reference/similarityMat.md)
[`similarityMatC`](http://sistm.github.io/NPflow/reference/similarityMatC.md)

## Author

Boris Hejblum
