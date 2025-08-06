
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `NPflow`

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/NPflow)](https://cran.r-project.org/package=NPflow)
[![R-CMD-check](https://github.com/sistm/NPflow/workflows/R-CMD-check/badge.svg)](https://github.com/sistm/NPflow/actions)
[![Downloads](https://cranlogs.r-pkg.org/badges/NPflow?color=blue)](https://www.r-pkg.org/pkg/NPflow)
<!-- badges: end -->

## Overview

`NPflow` is a package for performing **Bayesian estimation of Dirichlet
process mixtures of multivariate skew $t$-distributions**. It is
especially oriented towards flow-cytometry data preprocessing
applications, but can be useful for numerous other applications.

The main function of the package is `DPMpost()`.

The method implemented in this package is detailed in the following
article:

> Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019).
> Sequential Dirichlet Process Mixtures of Multivariate Skew
> t-distributions for Model-based Clustering of Flow Cytometry Data.
> *The Annals of Applied Statistics*, **13**(1):638-660. [\<doi:
> 10.1214/18-AOAS1209\>](https://doi.org/10.1214/18-AOAS1209) [\<arXiv:
> 1702.04407\>](https://arxiv.org/abs/1702.04407)

## Installation

The easiest way to get `NPflow` is to install it from
[CRAN](https://cran.r-project.org/package=NPflow):

``` r
install.packages("NPflow")
```

Or to get the development version from
[GitHub](https://github.com/sistm/NPflow):

``` r
#install.packages("devtools")
devtools::install_github("sistm/NPflow", ref="CRANrelease")
```

– Boris Hejblum
