# C++ implementation

C++ implementation

## Usage

``` r
similarityMat_nocostC(cc)
```

## Arguments

- cc:

  a matrix whose columns each represents a ()MCMC) partition

## Examples

``` r
c <- list(c(1,1,2,3,2,3), c(1,1,1,2,3,3),c(2,2,1,1,1,1))
similarityMat_nocostC(sapply(c, "["))
#> $similarity
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 1.0000000 1.0000000 0.3333333 0.0000000 0.0000000 0.0000000
#> [2,] 1.0000000 1.0000000 0.3333333 0.0000000 0.0000000 0.0000000
#> [3,] 0.3333333 0.3333333 1.0000000 0.3333333 0.6666667 0.3333333
#> [4,] 0.0000000 0.0000000 0.3333333 1.0000000 0.3333333 0.6666667
#> [5,] 0.0000000 0.0000000 0.6666667 0.3333333 1.0000000 0.6666667
#> [6,] 0.0000000 0.0000000 0.3333333 0.6666667 0.6666667 1.0000000
#> 

c2 <- list()
for(i in 1:10){
    c2 <- c(c2, list(rmultinom(n=1, size=1000, prob=rexp(n=1000))))
}

c3 <- sapply(c2, "[")

if(require(microbenchmark)){
library(microbenchmark)
microbenchmark(similarityMat(c3), similarityMat_nocostC(c3), times=2L)
}else{
cat("package 'microbenchmark' not available\n")
}
#> Unit: milliseconds
#>                       expr        min         lq       mean     median
#>          similarityMat(c3) 121.928105 121.928105 123.454694 123.454694
#>  similarityMat_nocostC(c3)   7.193267   7.193267   7.926916   7.926916
#>          uq        max neval
#>  124.981283 124.981283     2
#>    8.660565   8.660565     2
```
