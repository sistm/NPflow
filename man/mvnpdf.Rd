\name{mvnpdf}
\alias{mvnpdf}
\title{multivariate-Normal  probability density function}
\usage{
mvnpdf(x, mean, varcovM)
}
\arguments{
  \item{x}{data either a matrix}

  \item{mean}{mean vector or list of mean vectors (either a
  vector, a matrix or a list)}

  \item{varcovM}{}
}
\description{
multivariate-Normal probability density function
}
\examples{
mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
dnorm(1.96)

mvnpdf(x=matrix(rep(1.96,2), nrow=1, ncol=2),
      mean=c(0, 0), varcovM=diag(2)
)
}

