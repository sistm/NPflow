#'meanltivariate  Normal  probability density function
#'
#'@param x
#'
#'@param mean
#'
#'@param varcovM
#'
#'@export
#'
#'@example{
#'
#'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#'dnorm(1.96)
#'
#'mvnpdf(x=matrix(rep(1.96,2), nrow=1, ncol=2), 
#'       mean=c(0, 0), varcovM=diag(2)
#')
#'
#'}
#'

mvnpdf <- function(x, mean, varcovM){
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    if(is.null(mean)){
        mean <- rep(0, p)
        x0 <- x
    } else if(is.vector(mean) & length(mean)==p){
        x0 <- x-mean
    } else{
        stop("wrong input for mean")
    }
    
    
    if(is.null(varcovM)){
        varcovM <- diag(1, nrow=p, ncol=p)
    }
    if(dim(varcovM)[1]!=dim(varcovM)[2] | dim(varcovM)[1]!=p){
        stop("varcovM is not a square matrix or is of the wrong size")
    }
    
    R = chol(varcovM)
    xRinv <- x0%*%solve(R)
    logSqrtDetvarcovM <- sum(log(diag(R)))
    
    quadform <- xRinv%*%t(xRinv)
    y <- exp(-0.5*quadform - logSqrtDetvarcovM -p*log(2*pi)/2)
    
    return(y)
    
}