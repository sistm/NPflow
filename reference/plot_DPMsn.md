# Plot of a Dirichlet process mixture of skew normal distribution partition

Plot of a Dirichlet process mixture of skew normal distribution
partition

## Usage

``` r
plot_DPMsn(
  z,
  c,
  i = "",
  alpha = "?",
  U_SS,
  dims2plot = 1:nrow(z),
  ellipses = ifelse(length(dims2plot) < 3, TRUE, FALSE),
  gg.add = list(theme()),
  nbsim_dens = 1000
)
```

## Arguments

- z:

  data matrix `d x n` with `d` dimensions in rows and `n` observations
  in columns.

- c:

  allocation vector of length `n` indicating which observation belongs
  to which clusters.

- i:

  current MCMC iteration number.

- alpha:

  current value of the DP concentration parameter.

- U_SS:

  a list containing `"xi"`, `"psi"`, `"S"`, and `"df"`.

- dims2plot:

  index vector, subset of `1:d` indicating which dimensions should be
  drawn. Default is all of them.

- ellipses:

  a logical flag indicating whether ellipses should be drawn around
  clusters. Default is `TRUE` if only 2 dimensions are plotted, `FALSE`
  otherwise.

- gg.add:

  A list of instructions to add to the `ggplot2` instruction (see
  `gg-add`). Default is `list(theme())`, which adds nothing to the plot.

- nbsim_dens:

  number of simulated points used for computing clusters density
  contours in 2D plots. Default is `1000` points.

## Author

Boris Hejblum
