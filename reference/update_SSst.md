# Return updated sufficient statistics S for skew t-distribution with data matrix z

For internal use only.

## Usage

``` r
update_SSst(z, S, ltn, scale, df, hyperprior = NULL)
```

## Arguments

- z:

  data matrix

- S:

  previous sufficient statistics

- ltn:

  random effects

- df:

  skew t degrees of freedom

- hyperprior:

  Default is `NULL`
