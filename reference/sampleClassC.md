# C++ implementation of the multinomial sampling from a matrix of column vectors, each containing the sampling probabilities for their respective draw

C++ implementation of the multinomial sampling from a matrix of column
vectors, each containing the sampling probabilities for their respective
draw

## Usage

``` r
sampleClassC(probMat, Log = FALSE)
```

## Arguments

- probMat:

  a numeric matrix of dim `k x n` of containing column vectors of
  sampling probabilities for each class `k`.

- Log:

  a logical flag indicating whether the provided `probMat` is on the log
  scale or natural probability scale. Default is `FALSE` in which case
  it is considered on the natural probability scale.

## Value

a vector of integer of length `n` containing the multinomial draws for
each observation, i.e. the class allocation.
