# C++ implementation of the F-measure computation

C++ implementation of the F-measure computation

## Usage

``` r
FmeasureC(pred, ref)
```

## Arguments

- pred:

  vector of a predicted partition

- ref:

  vector of a reference partition

## Examples

``` r
pred <- c(1,1,2,3,2,3)
ref <- c(2,2,1,1,1,3)
FmeasureC(pred, ref)
#> [1] 0.8444444
```
