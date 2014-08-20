#'multivariate Skew-t  probability density function
#'
#'
#'@param x
#'
#'@param xi
#'
#'@param sigma
#'
#'@param psi
#'
#'@param df
#'
#'@export
#'
#'
mvstpdf <- function(x, xi, sigma, psi, df){
    
    if(is.null(x) | is.null(xi) | is.null(sigma) | is.null(psi)| is.null(df)){
        stop("some arguments are empty")
    }
    
    if(!is.matrix(x)){
        stop("x should be a matrix")
    }
    n <- dim(x)[2]
    p <- dim(x)[1]
    
    if(!is.list(xi)){
        if(is.null(xi)){
            stop("xi is empty")
        } else if(is.vector(xi) && length(xi)==p){
            x0 <- x-xi
        } else if(is.matrix(xi) && ncol(xi)==n){
            x0 <- x-xi
        } else{
            stop("wrong input for xi")
        }
    }else{
        x0 <- lapply(xi, function(v){x - v})
    }
    
    
    if(is.matrix(sigma)){
        #recovering original paremters
        omega <- sigma + tcrossprod(psi)
        omegaInv <- solve(omega)
        smallomega <- diag(sqrt(diag(omega)))
        alph <- (smallomega%*%omegaInv%*%psi
                 /as.vector(sqrt(1-crossprod(psi,omegaInv)%*%psi)))
        Qy <- apply(X=x0, MARGIN=2, FUN=function(v){crossprod(v,omegaInv)%*%v})
        
        if(dim(omega)[1]!=dim(omega)[2]){
            stop("omega is not a square matrix")
        }
        if(dim(omega)[1]!=p){
            stop("omega is of the wrong size")
        }        
        part1 <- 2*mvtpdf(x, mean=xi, varcovM=omega)
        part2 <- pt(q=(t(alph)%*%diag(1/sqrt(diag(omega)))%*%(x0)*
                        sqrt((df+p)/(df+Qy))),
                    df=df+p)
    }
    else{
        if(!is.list(sigma)){
            sigma <- apply(X=sigma, MARGIN=3, FUN=list)
            sigma <- lapply(sigma, FUN='[[', 1)
            x0 <- apply(X=x0, MARGIN=2, FUN=list)
            x0 <- lapply(x0, FUN='[[', 1)
            psi <- apply(X=psi, MARGIN=2, FUN=list)
            psi <- lapply(psi, FUN='[[', 1)
        }
        
        omega <- mapply(FUN=function(s,ps){s + tcrossprod(ps)},
                        s=sigma, ps=psi, SIMPLIFY=FALSE)
        omegaInv <- lapply(X=omega, FUN=solve)
        alph <- mapply(FUN=function(o, oI, ps){
            diag(sqrt(diag(o)))%*%oI%*%ps/sqrt(1-crossprod(ps,oI)%*%ps)[1,1]},
            o=omega, oI=omegaInv, ps=psi, SIMPLIFY=FALSE
        )
        Qy <- mapply(FUN=function(x, oI){
            apply(X=x, MARGIN=2,FUN=function(v){crossprod(v,oI)%*%v})},
            x=x0, oI=omegaInv, SIMPLIFY=FALSE)
        part1 <- 2*mvtpdf(x, mean=xi, varcovM=omega, df=df)
        part2 <- mapply(FUN=function(a, o, Q, d, x){
            pt(q=(crossprod(a,diag(1/sqrt(diag(o))))%*%(x)*
                sqrt((d+p)/(d+Q))),
        df=d+p)}, 
            x=x0, o=omega, Q=Qy, a=alph, d=df)
        
    }
    
    return(part1*part2)
    
}

