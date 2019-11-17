#' Sample from a Gamma distribution with a small shape parameter
#' 
#'When concentration parameter alpha is close to zero (only 1 non empty 
#'cluster being sampled), the data augmentation Gibbs sampling from Escobar 
#'& West relies on the sampling of from a Gamma distribution 
#'
#'For internal use.
#'
#'@param n the number of observations to sample
#'
#'@param shape the shape parameter of the gamma distribution. Default is 1.
#'
#'@param rate the shape parameter of the gamma distribution. Default is 1.
#'
#'@param do.log logical indicating whether the logarithm should be returned 
#'instead of the sampled value natural scale 
#'
#'@importFrom stats rgamma
#'
#'@keywords internal
#'
#'@reference Liu, C., Martin, R., & Syring, N. (2017). Efficient simulation 
#'from a gamma distribution with small shape parameter. 
#'Computational Statistics, 32(4), 1767–1775. 
#'https://doi.org/10.1007/s00180-016-0692-0
#'
#'
#'@source adapted from rgamss.R v2 (04/07/2015) by R. Martin 
#'https://www4.stat.ncsu.edu/~rmartin/Codes/rgamss.R 
#'
#'@examples
#'shape <- 0.001
#'log(rgamma(n = 10, shape = shape))
#'NPflow:::rgamma_ss(n = 10, shape = shape)
#'
rgamma_ss <- function(n = 1, shape = 1, rate = 1, do.log = TRUE){
  
  a <- shape
  if(a > 0.2){
    res_final <- stats::rgamma(n = n, shape = a, rate = rate)
  }else{
    #e1 <- 2.71828182845905 #sprintf("%.100f", exp(1))
    eta <- function(z, a, L){
      w <- a/(exp(1)*(1 - a))
      return(ifelse(z >= 0, exp(-z), w*L*exp(L*z)))
    }
    h <- function(z, a){
      return(exp(-z - exp(-z/a)))
    }
    rh <- function(a){
      L <- 1/a - 1
      ww <- 1/(1 + a/(exp(1)*(1 - a)))
      repeat{
        U <- runif(1)
        z <- ifelse(U <= ww, -log(U/ww), log(runif(1))/L)
        if(h(z, a)/eta(z, a, L) > runif(1)){ 
          return(z)
        }
      }
    }
    
    Z <- sapply(as.list(rep(a, n)), rh)
    res_temp <- - log(rate) - Z/a
    
    if(!do.log) {
      res_final <- exp(res_temp)
      if(any(res_final == 0)) {
        message("Insufficient numerical precision on the natural scale. Consider using the logarithm scale...")
      }
    }else{
      res_final <- res_temp
    }
  }
  return(res_final)
}


