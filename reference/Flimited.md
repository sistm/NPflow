# Compute a limited F-measure

A limited version of F-measure that only takes into account small
clusters

## Usage

``` r
Flimited(n_small_clst, pred, ref)
```

## Arguments

- n_small_clst:

  an integer for limit size of the small cluster

- pred:

  vector of a predicted partition

- ref:

  vector of a reference partition

## References

Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019)
Sequential Dirichlet Process Mixtures of Multivariate Skew
t-distributions for Model-based Clustering of Flow Cytometry Data. The
Annals of Applied Statistics, 13(1): 638-660. \<doi:
10.1214/18-AOAS1209\> \<arXiv: 1702.04407\>
<https://arxiv.org/abs/1702.04407>
[doi:10.1214/18-AOAS1209](https://doi.org/10.1214/18-AOAS1209)

## Examples

``` r
pred <- c(rep(1, 5),rep(2, 8),rep(3,10))
ref <- c(rep(1, 5),rep(c(2,3), 4),rep(c(3,2),5))
FmeasureC(pred, ref)
#> [1] 0.6292906
Flimited(6, pred, ref)
#> [1] 1
```
