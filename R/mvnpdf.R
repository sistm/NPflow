#'multivariate-Normal  probability density function
#'
#'@param x data either a matrix
#'
#'@param mean mean vector or list of mean vectors (either a vector, 
#'a matrix or a list)
#'
#'@param varcovM
#'
#'@export mvnpdf
#'
#'@examples
#'
#'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#'dnorm(1.96)
#'
#'mvnpdf(x=matrix(rep(1.96,2), nrow=1, ncol=2), 
#'       mean=c(0, 0), varcovM=diag(2)
#')
#'
#'
#'

mvnpdf <- function(x, mean, varcovM){
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    n <- dim(x)[2]
    p <- dim(x)[1]
   
    if(!is.list(mean)){
        if(is.null(mean)){
            stop("mean is empty")
        } else if(is.vector(mean) && length(mean)==p){
            x0 <- x-mean
        } else if(is.matrix(mean) && ncol(mean)==n){
            x0 <- x-mean
        } else{
            stop("wrong input for mean")
        }
    }else{
        if(!is.list(x)){
            x <- apply(X=x, MARGIN=2, FUN=list)
            x <- lapply(x, FUN='[[', 1)
        }
        x0 <- mapply('-', x=x, y=mean, SIMPLIFY=FALSE)
        
    }
    
    
    if(is.null(varcovM)){
        stop("varcovM is empty")
    }
    
    
    if(is.matrix(varcovM)){
        if(dim(varcovM)[1]!=dim(varcovM)[2]){
            stop("varcovM is not a square matrix")
        }
        if(dim(varcovM)[1]!=p){
            stop("varcovM is of the wrong size")
        }
        
        if(is.vector(x0)){
            x0 <- matrix(x0, ncol=1)
        }
        
        Rinv = backsolve(chol(varcovM),diag(p))
        xRinv <- apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv)
        logSqrtDetvarcovM <- sum(log(diag(Rinv)))
        
        quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
        y <- exp(-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)
        
#         dMvn <- function(X,mu,Sigma) {
#             k <- ncol(X)
#             rooti <- backsolve(chol(Sigma),diag(k))
#             quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
#             return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
#         }
        
    }else{
        if(!is.list(varcovM)){
            varcovM <- apply(X=varcovM, MARGIN=3, FUN=list)
            varcovM <- lapply(varcovM, FUN='[[', 1)
        }    
        if(!is.list(x0)){
            x0 <- apply(X=x0, MARGIN=2, FUN=list)
            x0 <- lapply(x0, FUN='[[', 1)
        }
        R <- lapply(varcovM, FUN=chol)
        Rinv <- lapply(X=R, FUN=solve)
        xRinv <- mapply(FUN='%*%', x=x0, y=Rinv, SIMPLIFY=FALSE)
        logSqrtDetvarcovM <- lapply(X=R, FUN=function(X){sum(log(diag(X)))})
        quadform <- lapply(X=xRinv, FUN=tcrossprod)
        y <- mapply(FUN=function(x,y){exp(-0.5*x - y -p*log(2*pi)/2)},
                    x=quadform, y=logSqrtDetvarcovM)  
    }
    
    
    return(y)
    
}