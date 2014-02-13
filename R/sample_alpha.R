sample_alpha <- function(alpha_old, n, K, a, b){
    
    # Sample scale factor in Dirichlet Process 
    # Cf data augmentation by West et al., 1992
    # n is the number of data points 
    # K the number of clusters
    x <- rbeta(n=1, shape1=alpha_old + 1, shape2=n)
    temp <- (a+K-1) / (n*(b-log(x)))
    pi = temp/(1+temp);
    u=runif(1)
    if (u<pi){
        gamma_new = rgamma(n=1, shape=a + K, scale=1/(b - log(x)))
    } else{
        gamma_new = rgamma(n=1, shape=a + K - 1, scale=1/(b - log(x)))
    }
    return(gamma_new)
}