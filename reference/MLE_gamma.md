# MLE for Gamma distribution

Maximum likelihood estimation of Gamma distributed observations
distribution parameters

## Usage

``` r
MLE_gamma(g)
```

## Arguments

- g:

  a list of Gamma distributed observation.

## Examples

``` r

g_list <- list()
for(i in 1:1000){
 g_list <- c(g_list, rgamma(1, shape=100, rate=5))
}

mle <- MLE_gamma(g_list)
#> Warning: Unable to estimate Gamma hyperpriors properly
#>             (this can happen when only a few clusters are fitted).
#>             => Non informative values are returned instead
mle
#> $shape
#> [1] 1e-04
#> 
#> $scale
#> [1] 1e-04
#> 
#> $rate
#> [1] 10000
#> 
```
