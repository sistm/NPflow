# Convergence diagnostic plots

Convergence diagnostic plots

## Usage

``` r
plot_ConvDPM(
  MCMCsample,
  from = 1,
  to = length(MCMCsample$logposterior_list),
  shift = 0,
  thin = 1,
  ...
)
```

## Arguments

- MCMCsample:

  a `DPMMclust` or `summaryDPMMclust` object.

- from:

  the MCMC iteration from which the plot should start. Default is `1`.

- to:

  the MCMC iteration up until which the plot should stop. Default is
  `1`.

- shift:

  a number of initial iterations not to be displayed. Default is `0`.

- thin:

  integer giving the spacing at which MCMC iterations are kept. Default
  is `1`, i.e. no thining.

- ...:

  further arguments passed to or from other methods
