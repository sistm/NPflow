# Sample from a Wishart distribution
# n : degrees of freedom
# lambda : scale parameter

wishrnd <- function(n, lambda){
  
  p <- dim(lambda)[1]
  p2 <- dim(lambda)[2]
  
  if(p!=p2){
    stop('Error : Matrix not square\n')
  }
  
  x=chol(lambda)%*%matrix(rnorm(p*n), nrow=p, ncol=n)
  W=x%*%t(x)
  
  return(W)
}