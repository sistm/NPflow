# Sample from an inverse Wishart distribution
# n : degrees of freedom
# lambda : scale parameter

invwishrnd <- function(n,lambda){
  iS=wishrnd(n = n, lambda = solve(lambda))
  S=solve(iS)
  return(S)
}