# Methods for a summary of a `DPMMclust` object

Methods for a summary of a `DPMMclust` object

## Usage

``` r
# S3 method for class 'summaryDPMMclust'
print(x, ...)

# S3 method for class 'summaryDPMMclust'
plot(
  x,
  hm = FALSE,
  nbsim_densities = 5000,
  hm_subsample = NULL,
  hm_order_by_clust = TRUE,
  gg.add = list(theme_bw()),
  ...
)
```

## Arguments

- x:

  a `summaryDPMMclust` object.

- ...:

  further arguments passed to or from other methods

- hm:

  logical flag to plot the heatmap of the similarity matrix. Default is
  `FALSE`.

- nbsim_densities:

  the number of simulated observations to be used to plot the density
  lines of the clusters in the point estimate partition plot

- hm_subsample:

  a integer designating the number of observations to use when plotting
  the heatmap. Used only if `hm` is `TRUE`. \#'Default is `NULL` in
  which no subsampling is done and all observations are plotted.

- hm_order_by_clust:

  logical flag indicating whether observations should be ordered
  according to the point estimate first. Used only if `hm` is `TRUE`.
  Default is `TRUE`.

- gg.add:

  a list of instructions to add to the `ggplot2` instruction (see
  `gg-add`). Default is `list(theme())`, which adds nothing to the plot.

## Author

Boris Hejblum
