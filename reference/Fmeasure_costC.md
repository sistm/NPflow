# Multiple cost computations with the F-measure as the loss function

C++ implementation of multiple cost computations with the F-measure as
the loss function using the Armadillo library

## Usage

``` r
Fmeasure_costC(c)
```

## Arguments

- c:

  a matrix where each column is one MCMC partition

## Value

a list with the following elements:

- `Fmeas:`:

  TODO

- `cost:`:

  TODO

## Examples

``` r
library(NPflow)
c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
#Fmeasure_costC(sapply(c, "["))

if(interactive()){
c2 <- list()
for(i in 1:100){
    c2 <- c(c2, list(rmultinom(n=1, size=2000, prob=rexp(n=2000))))
}
Fmeasure_costC(sapply(c2, "["))
}
```
