#'Bayesian Nonparametrics for Automatic Gating of Flow Cytometry data
#'
#'Dirichlet process mixture of multivariate normal, skew normal or skew t-distributions
#'modeling oriented towards flow-cytometry data pre-processing applications.
#'
#'
#'\tabular{ll}{
#'Package: \tab NPflow\cr
#'Type: \tab Package\cr
#'Version: \tab 0.13.5\cr
#'Date: \tab 2024-01-13\cr
#'License:\tab \href{http://www.gnu.org/licenses/lgpl.txt}{LGPL-3}\cr
#'}
#'The main function in this package is \code{\link{DPMpost}}.
#'
#'@author Boris P. Hejblum, Chariff Alkhassim, Francois Caron
#'--- Maintainer: Boris P. Hejblum
#'
#'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F and Thiebaut R (2019) 
#'Sequential Dirichlet Process Mixtures of Multivariate Skew t-distributions for 
#'Model-based Clustering of Flow Cytometry Data. The Annals of Applied Statistics, 
#'13(1): 638-660. <doi: 10.1214/18-AOAS1209> <arXiv: 1702.04407> 
#'\url{https://arxiv.org/abs/1702.04407} \doi{10.1214/18-AOAS1209}
#'
#'@name NPflow-package
#'@aliases NPflow
#'
#'@useDynLib NPflow, .registration = TRUE
#'@importFrom Rcpp evalCpp
#'
"_PACKAGE"