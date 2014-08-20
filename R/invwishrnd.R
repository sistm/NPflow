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
invwishrnd <- function(n,lambda){
    iS=wishrnd(n = n, Sigma = solve(lambda))
    S=solve(iS)
    return(S)
}