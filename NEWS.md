# News about the `NPflow` R package

### Main changes in Version 0.13.6 (2025-08-20) --- *this is only a minor release*:
* added a pkgdown website
* added `MASS` in Suggest as it is no longer imported by `ggplot2`

### Main changes in Version 0.13.5 (2024-01-13) --- *this is only a minor release*:
* bug fixed in try-catch syntax

### Main changes in Version 0.13.3 (2020-02-06) --- *this is only a minor release*:
* bug fixed in documentation and compiling standards

### Main changes in Version 0.13.1 (2017-08-02) --- *this is only a minor release*:
* updated documentation and compiling standards


### Main changes in Version 0.13.0 (2017-04-10):
* use of `itertools` package for parallel computations
* use of log-scale probabilities in slice-samplers for improved precision in computation
* *bug fix*: parallel settings did not work in `DPMpost` (thanks to **Mike Jiang** for noticing this)
* updated documentation


### Main changes in Version 0.12.0 (2017-04-04) --- *this is only a minor release*:
* registration of compiled functions
* updated help and citation info


### Main changes in Version 0.10.0 (2016-05-02):
* optimization of the log-posterior computation for normal and skew normal distributions
* improvement of the parallel sampler for normal distributions
* added arguments in `DPMGibs...` functions to allow more informative priors if needed
* *bug fix*: log of integers disambiguated in C++ function (compilation error on Solaris)

