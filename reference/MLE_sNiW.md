# MLE for sNiW distributed observations

Maximum likelihood estimation of Normal inverse Wishart distributed
observations

## Usage

``` r
MLE_sNiW(xi_list, psi_list, S_list, doPlot = TRUE)
```

## Arguments

- xi_list:

  a list of length `N` whose elements are observed vectors of length `d`
  of the mean parameters xi.

- psi_list:

  a list of length `N` whose elements are observed vectors of length `d`
  of the skew parameters psi.

- S_list:

  a list of length `N` whose elements are observed variance-covariance
  matrices of dimension `d x d`.

- doPlot:

  a logical flag indicating whether the algorithm progression should be
  plotted. Default is `TRUE`.

## Author

Boris Hejblum, Chariff Alkhassim

## Examples

``` r
hyperG0 <- list()
hyperG0$b_xi <- c(0.3, -1.5)
hyperG0$b_psi <- c(0, 0)
hyperG0$kappa <- 0.001
hyperG0$D_xi <- 100
hyperG0$D_psi <- 100
hyperG0$nu <- 35
hyperG0$lambda <- diag(c(0.25,0.35))

xi_list <- list()
psi_list <- list()
S_list <- list()
for(k in 1:1000){
 NNiW <- rNNiW(hyperG0, diagVar=FALSE)
 xi_list[[k]] <- NNiW[["xi"]]
 psi_list[[k]] <- NNiW[["psi"]]
 S_list[[k]] <- NNiW[["S"]]
}

mle <- MLE_sNiW(xi_list, psi_list, S_list)
mle
#> $U_xi
#> [1]  0.268813 -1.496686
#> 
#> $U_psi
#> [1] -0.036632007  0.003716794
#> 
#> $U_B
#>               [,1]          [,2]
#> [1,]  9.952192e-03 -7.728564e-05
#> [2,] -7.728564e-05  9.834377e-03
#> 
#> $U_df
#> [1] 34.82543
#> 
#> $U_Sigma
#>              [,1]         [,2]
#> [1,] 0.2448164671 0.0006174342
#> [2,] 0.0006174342 0.3532194135
#> 
```
