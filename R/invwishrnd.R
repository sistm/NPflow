#' Sample from a inverse-Wishart distribution
#' 
#' For internal use only.
#' 
#'@param n degrees of freedom
#'
#'@param lambda scale parameter
#' 
#'@keywords internal
#'
#'@export
#'
#'
invwishrnd <- function(n, lambda, lambda_solved){
  S <- solve(wishrnd(n = n, Sigma = lambda_solved))
  return(S)
}
