#'multivariate Student's t-distribution probability density function
#'
#'
#'@param x data matrix
#'
#'@param mean mean vector or list of mean vectors (either a vector, 
#'a matrix or a list)
#'
#'@param varcovM
#'
#'@param df
#'
#'@export
#'
#'@examples
#'mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10000000)
#'mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1))
#'
#'mvtpdf(x=matrix(1.96), mean=0, varcovM=diag(1), df=10)
#'
#'mvtpdf(x=matrix(rep(1.96,2), nrow=2, ncol=1), 
#'       mean=c(0, 0), varcovM=diag(2), df=10
#')
#'
mvtpdf <- function(x, mean, varcovM, df){
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    n <- ncol(x)
    p <- nrow(x)
    
    
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
        
        x0 <- lapply(mean, function(v){x - v})
        
    }
    
    if(is.null(varcovM)){
        stop("varcovM is empty")
    }


    if(is.matrix(varcovM)){
        if(ncol(varcovM)!=nrow(varcovM)){
            stop("varcovM is not a square matrix")
        }
        if(nrow(varcovM)!=p){
            stop("varcovM is of the wrong size")
        }
        
        if(is.vector(x0)){
            x0 <- matrix(x0, ncol=1)
        }
        
       
        
        
        Rinv = backsolve(chol(varcovM),x=diag(p))
        xRinv <- matrix(apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv))
        logSqrtDetvarcovM <- sum(log(diag(Rinv)))
        a <- lgamma((df + p)/2)-lgamma(df/2)-p/2*log(df*pi)
        quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
        y <- (1+quadform/df)^(-(df+p)/2)*exp(a+logSqrtDetvarcovM)
        
#         dMvn <- function(X,mu,Sigma) {
#             k <- ncol(X)
#             rooti <- backsolve(chol(Sigma),diag(k))
#             quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
#             return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
#         }
    }else if(!is.list(mean)){
        if(!is.list(varcovM)){
            varcovM <- apply(X=varcovM, MARGIN=3, FUN=list)
            varcovM <- lapply(varcovM, FUN='[[', 1)
        }    
        if(!is.list(x0)){
            x0 <- apply(X=x0, MARGIN=2, FUN=list)
            x0 <- lapply(x0, FUN='[[', 1)
        }    
        if(!is.list(df)){
            df <- lapply(df, FUN='[')
        }
        
        likelihood <- function(x0, varcovM, df){
            p <- length(x0)            
            Rinv = backsolve(chol(varcovM),x=diag(p))
            xRinv <- x0 %*% Rinv
            logSqrtDetvarcovM <- sum(log(diag(Rinv)))
            a <- lgamma((df + p)/2)-lgamma(df/2)-p/2*log(df*pi)
            quadform <- tcrossprod(xRinv)
            y <- (1+quadform/df)^(-(df+p)/2)*exp(a+logSqrtDetvarcovM)
        }
        y <-mapply(FUN=likelihood, x0, varcovM, df)
        
        
     }else{
        R <- lapply(varcovM, FUN=chol)
        Rinv <- lapply(X=R, FUN=backsolve, x=diag(p))
        xRinv <- mapply(FUN=function(x,y){apply(X=x,MARGIN=2, FUN=crossprod, y=y)}, x=x0, y=Rinv, SIMPLIFY=FALSE)
        logSqrtDetvarcovM <- lapply(X=Rinv, FUN=function(X){sum(log(diag(X)))})
        a <- lapply(df, FUN=function(x){lgamma((x + p)/2)-lgamma(x/2)-p/2*log(x*pi)})
        quadform <- lapply(X=xRinv, FUN=function(x){apply(X=x, MARGIN=2, FUN=crossprod)})
        y <- mapply(FUN=function(u,v,w,z){(1+v/u)^(-(u+p)/2)*exp(w+z)},
                    u=df, v=quadform, w=a, z=logSqrtDetvarcovM)  
    }
    
    return(y)
    
}