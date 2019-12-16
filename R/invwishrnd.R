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
invwishrnd <- function(n, lambda_solved){
  p <- ncol(lambda_solved)
  S <- try(solve(wishrnd(n = n, Sigma = lambda_solved)), silent=TRUE)
  if(inherits(S, "try-error")){
    S <- solve(wishrnd(n = n, Sigma = solve((lambda + diag(p)))))
  }
  return(S)
}
