# C++ implementation of similarity matrix computation using pre-computed distances

C++ implementation of similarity matrix computation using pre-computed
distances

## Usage

``` r
NuMatParC(c, d)
```

## Arguments

- c:

  an MCMC partitions of length `n`.

- d:

  a symmetric `n x n` matrix containing distances between each group
  distributions.

## Author

Boris Hejblum, Chariff Alkhassim

## Examples

``` r
c <- c(1,1,2,3,2,3)
d <- matrix(runif(length(c)^2),length(c))
NuMatParC(c,d)
#> $NuMatParC
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 0.0000000 0.0000000 0.9259166 0.2847929 0.9259166 0.2847929
#> [2,] 0.0000000 0.0000000 0.9259166 0.2847929 0.9259166 0.2847929
#> [3,] 0.9259166 0.9259166 0.0000000 0.2237126 0.0000000 0.2237126
#> [4,] 0.2847929 0.2847929 0.2237126 0.0000000 0.8303850 0.0000000
#> [5,] 0.9259166 0.9259166 0.0000000 0.8303850 0.0000000 0.2237126
#> [6,] 0.2847929 0.2847929 0.2237126 0.0000000 0.2237126 0.0000000
#> 

```
