# Sample from a Wishart distribution
# n : degrees of freedom
# Sigma : scale parameter

wishrnd <- function(n, Sigma){
  
  p <- dim(Sigma)[1]
  p2 <- dim(Sigma)[2]
  
  if(p!=p2){
    stop('Error : Matrix not square\n')
  }
  
  x=chol(Sigma)%*%matrix(rnorm(p*n), nrow=p, ncol=n)
  W=x%*%t(x)
  
  return(W)
}