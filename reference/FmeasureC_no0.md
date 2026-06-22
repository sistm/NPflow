# C++ implementation of the F-measure computation without the reference class 0

Aghaeepour in FlowCAP 1 ignore the reference class labeled "0"

## Usage

``` r
FmeasureC_no0(pred, ref)
```

## Arguments

- pred:

  vector of a predicted partition

- ref:

  vector of a reference partition

## References

N Aghaeepour, G Finak, H Hoos, TR Mosmann, RR Brinkman, R Gottardo, RH
Scheuermann, Critical assessment of automated flow cytometry data
analysis techniques, *Nature Methods*, 10(3):228-38, 2013.

## Examples

``` r
library(NPflow)
pred <- c(1,1,2,3,2,3)
ref <- c(2,2,0,0,0,3)
FmeasureC(pred, ref)
#> [1] 0.8444444
FmeasureC_no0(pred, ref)
#> [1] 1
```
