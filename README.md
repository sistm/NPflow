
<!-- README.md is generated from README.Rmd. Please edit that file -->
`NPflow`
========

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NPflow)](https://cran.r-project.org/package=NPflow) [![Travis-CI Build Status](https://travis-ci.org/borishejblum/NPflow.svg?branch=CRANrelease)](https://travis-ci.org/borishejblum/NPflow) [![Downloads](https://cranlogs.r-pkg.org/badges/NPflow?color=blue)](https://www.r-pkg.org/pkg/NPflow)

Overview
--------

`NPflow` is a package for performing **Bayesian estimation of Dirichlet process mixtures of multivariate skew *t*-distributions**. It is especially oriented towards flow-cytometry data preprocessing applications, but can be useful for numerous other applications.

The main function of the package is `DPMpost()`.

The method implemented in this package is detailed in the following article:

> Hejblum BP, Alkhassim C, Gottardo R, Caron F, Thiebaut R, Sequential Dirichlet Process Mixtures of Multivariate Skew *t*-distributions for Model-based Clustering of Flow Cytometry Data, 2017, *submitted*. [arXiv:1702.04407](https://arxiv.org/abs/1702.04407v2).

Installation
------------

The easiest way to get `NPflow` is to install it from [CRAN](https://cran.r-project.org/package=NPflow):

``` r
install.packages("NPflow")
```

Or to get the development version from [GitHub](https://github.com/borishejblum/NPflow):

``` r
#install.packages("devtools")
devtools::install_github("borishejblum/NPflow")
```

-- Boris Hejblum
