#' Sample from a Wishart distribution
#'
#' For internal use only.
#'
#'@param n degrees of freedom
#'
#'@param Sigma scale parameter (matrix)
#'
#'@keywords internal
#'
#'@importFrom stats rnorm
#'
#'@export wishrnd

wishrnd <- function(n, Sigma){

  p <- nrow(Sigma)
  p2 <- ncol(Sigma)
  if(p!=p2){
    stop('scale parameter Sigma is not a square matrix')
  }
  
  x <- matrix(rnorm(n = n*p), nrow=n, ncol=p) %*% chol(Sigma)
  W <- crossprod(x)

  return(W)
}
