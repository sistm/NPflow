#'Sampler for the concentration parameter of a Dirichlet process
#'  
#'A Gamma prior is used
#'
#'@references West, M. (1992). Hyperparameter estimation in Dirichlet process mixture models. Technical Report, Duke University.
#'
#'@export
#'
#'@examples
#' #Test with a fixed K
#' ####################
#' 
#' alpha_init <- 1000
#' N <- 100000
#' #n=500
#' n=100000
#' K <- 80
#' a <- 0.0001
#' b <- a
#' alphas <- numeric(N)
#' alphas[1] <- alpha_init
#' for (i in 2:N){
#'  alphas[i] <- sample_alpha(alpha_old = alphas[i-1], n=n, K=K, a=a, b=b)
#' }
#' 
#' postalphas <- alphas[floor(N/2):N]
#' alphaMMSE <- mean(postalphas)
#' alphaMAP <- density(postalphas)$x[which.max(density(postalphas)$y)]
#' 
#' expK <- sum(alphaMMSE/(alphaMMSE+0:(n-1)))
#' round(expK)
#' 
#' 
#'  prioralpha <- data.frame("alpha"=rgamma(n=5000, a,1/b), 
#'                          "distribution" =factor(rep("prior",5000), 
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(prioralpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="red")
#'        + ggtitle(paste("Prior distribution on alpha: Gamma(", a, 
#'                  ",", b, ")\n", sep=""))
#'       )
#'  p
#' 
#' postalpha.df <- data.frame("alpha"=postalphas, 
#'                          "distribution" = factor(rep("posterior",length(postalphas)), 
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(postalpha.df, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..), binwidth=.1,
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of alpha\n")
#'        + geom_vline(aes(xintercept=mean(alpha, na.rm=T)),   # Ignore NA values for mean
#'                     color="red", linetype="dashed", size=1)  # Overlay with transparent density plot            
#'      )
#'  p
#' 
#'
#'
#'
#'
#'

sample_alpha <- function(alpha_old, n, K, a, b){
    
    # Sample scale factor in Dirichlet Process 
    # Cf data augmentation by West et al., 1992
    # n is the number of data points 
    # K the number of clusters
    x <- rbeta(n=1, shape1=alpha_old + 1, shape2=n)
    temp <- (a+K-1) / (n*(b-log(x)))
    pi <- temp/(1+temp)
    u <- runif(1)
    if (u<pi){
        gamma_new = rgamma(n=1, shape=a + K, scale=1/(b - log(x)))
    } else{
        gamma_new = rgamma(n=1, shape=a + K - 1, scale=1/(b - log(x)))
    }
    return(gamma_new)
}

