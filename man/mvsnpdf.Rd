\name{mvsnpdf}
\alias{mvsnpdf}
\title{multivariate Skew-Normal  probability density function}
\usage{
mvsnpdf(x, xi, sigma, psi)
}
\arguments{
  \item{x}{}

  \item{mean}{}

  \item{omega}{}
}
\description{
multivariate Skew-Normal probability density function
}
\examples{
mvsnpdf(x=matrix(1.96), mean=0, omega=diag(1))
dnorm(1.96)

mvsnpdf(x=matrix(rep(1.96,2), nrow=1, ncol=2),
      mean=c(0, 0), omega=diag(2)
)
}

