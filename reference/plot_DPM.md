# Plot of a Dirichlet process mixture of gaussian distribution partition

Plot of a Dirichlet process mixture of gaussian distribution partition

## Usage

``` r
plot_DPM(
  z,
  U_mu = NULL,
  U_Sigma = NULL,
  m,
  c,
  i,
  alpha = "?",
  U_SS = NULL,
  dims2plot = 1:nrow(z),
  ellipses = ifelse(length(dims2plot) < 3, TRUE, FALSE),
  gg.add = list(theme())
)
```

## Arguments

- z:

  data matrix `d x n` with `d` dimensions in rows and `n` observations
  in columns.

- U_mu:

  either a list or a matrix containing the current estimates of mean
  vectors of length `d` for each cluster. Default is `NULL` in which
  case `U_SS` has to be provided.

- U_Sigma:

  either a list or an array containing the `d x d` current estimates for
  covariance matrix of each cluster. Default is `NULL` in which case
  `U_SS` has to be provided.

- m:

  vector of length `n` containing the number of observations currently
  assigned to each clusters.

- c:

  allocation vector of length `n` indicating which observation belongs
  to which clusters.

- i:

  current MCMC iteration number.

- alpha:

  current value of the DP concentration parameter.

- U_SS:

  a list containing `"mu"` and `"S"`. Default is `NULL` in which case
  `U_mu` and `U_Sigma` have to be provided.

- dims2plot:

  index vector, subset of `1:d` indicating which dimensions should be
  drawn. Default is all of them.

- ellipses:

  a logical flag indicating whether ellipses should be drawn around
  clusters. Default is `TRUE` if only 2 dimensions are plotted, `FALSE`
  otherwise.

- gg.add:

  a list of instructions to add to the `ggplot2` instruction (see
  `gg-add`). Default is `list(theme())`, which adds nothing to the plot.

## Author

Boris Hejblum
