#' Sample from a normal inverse Wishart distribution 
#' whose parameter are given by the structure SufStat
#'
#' 
#' For internal use only.
#' 
#'@keywords internal
#'
#'@export nniw_rnd
#'

nniw_rnd <- function(SufStat){
  b0_xi = SufStat[["b_xi"]]
  b0_psi = SufStat[["b_psi"]]
  B = SufStat[["B"]] # B
  B0 <- diag(c(1/SufStat[["D_xi"]], 1/SufStat[["D_psi"]]))
  nu0 = SufStat[["nu"]] #c
  lambda0 = SufStat[["lambda"]] #C
  
  if(is.null(B)){
      B <- B0
  }
  
  # Sample S from an inverse Wishart distribution
  S = invwishrnd(n = nu0, lambda = lambda0)
  
  # Sample mu from a normal distribution
  d <- length(b0_xi)
  muSupp <- (c(b0_xi, b0_psi) 
             + matrix(rnorm(2*d), nrow=1, ncol=2*d)%*%chol(B%x%S))
  xi <- muSupp[1:d]
  psi <- muSupp[(d+1):(2*d)]
  
  
  return(list("S"=S, "xi"=xi, "psi"=psi))
}
