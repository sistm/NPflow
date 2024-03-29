% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postProcess.DPMMclust.R
\name{postProcess.DPMMclust}
\alias{postProcess.DPMMclust}
\title{Post-processing Dirichlet Process Mixture Models results to get
a mixture distribution of the posterior locations}
\usage{
postProcess.DPMMclust(
  x,
  burnin = 0,
  thin = 1,
  gs = NULL,
  lossFn = "F-measure",
  K = 10,
  ...
)
}
\arguments{
\item{x}{a \code{DPMMclust} object.}

\item{burnin}{integer giving the number of MCMC iterations to burn (defaults is half)}

\item{thin}{integer giving the spacing at which MCMC iterations are kept.
Default is \code{1}, i.e. no thining.}

\item{gs}{optional vector of length \code{n} containing the gold standard
partition of the \code{n} observations to compare to the point estimate.}

\item{lossFn}{character string specifying the loss function to be used.
Either "F-measure" or "Binder" (see Details). Default is "F-measure".}

\item{K}{integer giving the number of mixture components. Default is \code{10}.}

\item{...}{further arguments passed to or from other methods}
}
\value{
a \code{list}:
     \item{\code{burnin}:}{an integer passing along the \code{burnin} argument}
     \item{\code{thin}:}{an integer passing along the \code{thin} argument}
     \item{\code{lossFn}:}{a character string passing along the \code{lossFn} argument}
     \item{\code{point_estim}:}{}
     \item{\code{loss}:}{}
     \item{\code{index_estim}:}{}
}
\description{
Post-processing Dirichlet Process Mixture Models results to get
a mixture distribution of the posterior locations
}
\details{
The cost of a point estimate partition is calculated using either a pairwise
coincidence loss function (Binder), or 1-Fmeasure (F-measure).
}
\seealso{
\code{\link{similarityMat}} \code{\link{summary.DPMMclust}}
}
\author{
Boris Hejblum
}
