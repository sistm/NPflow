# C++ implementation

C++ implementation

## Usage

``` r
vclust2mcoclustC(c)
```

## Arguments

- c:

  is an MCMC partition

## Author

Chariff Alkhassim

## Examples

``` r
cc <- c(1,1,2,3,2,3)
vclust2mcoclustC(cc)
#> $Coclust
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    1    0    0    0    0
#> [2,]    1    0    0    0    0    0
#> [3,]    0    0    0    0    1    0
#> [4,]    0    0    0    0    0    1
#> [5,]    0    0    1    0    0    0
#> [6,]    0    0    0    1    0    0
#> 

```
