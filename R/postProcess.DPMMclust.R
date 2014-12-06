#'Post-processing Dirichlet Process Mixture Models results to get 
#'a mixture distribution of the posterior locations
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin integer giving the number of MCMC iterations to burn (defaults is half)
#'
#'@param thin integer giving the spacing at which MCMC iterations are kept. 
#'Default is \code{1}, i.e. no thining.
#'
#'@param lossFn character string specifying the loss function to be used.
#'Either "F-measure" or "Binder" (see Details). Default is "F-measure".
#'
#'@param gs optionnal vector of length \code{n} containing the gold standard 
#'partition of the \code{n} observations to compare to the point estimate
#'
#'@param ... further arguments passed to or from other methods
#'
#'@return a \code{list}: 
#'  \itemize{
#'      \item{\code{burnin}:}{an integer passing along the \code{burnin} argument}
#'      \item{\code{thin}:}{an integer passing along the \code{thin} argument}
#'      \item{\code{lossFn}:}{a character string passing along the \code{lossFn} argument}
#'      \item{\code{point_estim}:}{}
#'      \item{\code{loss}:}{}
#'      \item{\code{index_estim}:}{}
#'  }
#'
#'@details The cost of a point estimate partition is calculated using either a pairwise
#' coincidence loss function (Binder), or 1-Fmeasure (F-measure).
#'
#'@author Boris Hejblum
#'
#'@export 
#'
#'@importFrom gplots heatmap.2
#'
#'@seealso \link{similarityMat, summary.DPMMclust}
#'
postProcess.DPMMclust <- function(x, burnin=0, thin=1, gs=NULL, lossFn="F-measure", K=10,...){
    
    x_invar <- burn.DPMMclust(x, burnin = burnin, thin=thin)
    
    
    xi_list <- list()
    psi_list <- list()
    
    S_list <- list()
    
    #m_final <- list()
    #S_final <- list()
    
    for(i in 1:length(x_invar$U_SS_list)){
        xi_list <- c(xi_list, sapply(x_invar$U_SS_list[[i]], "[", "xi"))
        psi_list <- c(psi_list, sapply(x_invar$U_SS_list[[i]], "[", "psi"))
        
        S_list <- c(S_list, sapply(x_invar$U_SS_list[[i]], "[", "S"))
        
        #m_final <- c(m_final, 
        #             mapply(FUN=function(v1,v2){c(v1, v2)}, v1=xi_list, 
        #                    v2=psi_list, SIMPLIFY = FALSE)
        #)
        #B_list <- sapply(x_invar$U_SS_list[[i]], "[", "B")
        #S_final <- c(S_final, 
        #             mapply(FUN=function(M1,M2){M1%x%M2}, M1=B_list, M2=S_list, SIMPLIFY = FALSE)
        #)
        
    }
    
    if(x$clust_distrib=="skewT"){
        if(K>1){
            mle <- MLE_skewT_mmEM(xi_list, psi_list, S_list, 
                                  hyperG0 = x$hyperG0, K=K)
        }
        else{
            mle <- MLE_skewT(xi_list, psi_list, S_list)
        }
    }
    else{
        stop("clust_distrib is not skewT\n other distrib nort implemented yet")
    }
    
    return(mle)
}

#'EM MLE for mixture of sNiW
#'
#'Maximum likelihood estimation of mixture of 
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@rdname MLE_skewT_mmEM
#'
#'@export MLE_skewT_mmEM
#'
#'@examples
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(0.3, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 100
#'hyperG0$D_psi <- 100
#'hyperG0$nu <- 20
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'hyperG02 <- list()
#'hyperG02$b_xi <- c(-1, 2)
#'hyperG02$b_psi <- c(-0.1, 0.5)
#'hyperG02$kappa <- 0.001
#'hyperG02$D_xi <- 10
#'hyperG02$D_psi <- 10
#'hyperG02$nu <- 3
#'hyperG02$lambda <- 0.5*diag(2)
#'
#'for(k in 1001:2000){
#'  NNiW <- rNNiW(hyperG02, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'mle <- MLE_skewT_mmEM(xi_list, psi_list, S_list, hyperG0, K=2)
#'mle
#'
MLE_skewT_mmEM <- function( xi_list, psi_list, S_list, hyperG0, K, maxit=50, tol=1E-1){
    
    
    N <- length(xi_list)
    d <- length(hyperG0[[1]])
    
    if(length(psi_list) != N | length(S_list) != N){
        stop("Number of MCMC iterations not matching")
    }
    
    U_xi <- list() #matrix(0, nrow=d,ncol=K)
    U_psi <- list() #matrix(0, nrow=d,ncol=K)
    U_Sigma <- list() # array(dim=c(d,d,K))
    U_B <- list() #array(dim=c(2,2,K))
    U_df <- list() #numeric(K)
    
    
    #initialisation
    weights <- rep(1/K, K)
    for(k in 1:K){
        #sampling the cluster parameters
        NNiW <- rNNiW(hyperG0, diagVar=FALSE)
        U_xi[[k]] <- NNiW[["xi"]]
        U_psi[[k]] <- NNiW[["psi"]]
        U_Sigma[[k]] <- NNiW[["S"]]
        U_B[[k]] <- diag(0.01, 2)
        U_df[[k]] <- d+1
    }
    
    loglik <- numeric(maxit+1)
    loglik[1] <- -Inf
    Q <- numeric(maxit+1)
    Q[1] <- -Inf
    
    for(i in 1:maxit){
        
        r <- mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                          U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                          U_Sigma0 = U_Sigma, U_df0 = U_df)
        r <- apply(X=r, MARGIN=2, FUN=function(x){x+log(weights)})
        r <- apply(X=r, MARGIN=2, FUN=function(x){x - log(sum(exp(x)))})
        r[which(is.infinite(r))] <- -Inf
        r <- exp(r)
        
        
        
        #M step
        N_k <- rowSums(r)
        weights  <- N_k/N
        cat("weights:", weights, "\n")
        
        for(k in 1:K){
            U_xi[[k]] <- colSums(apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))/N_k[k]
            U_psi[[k]] <- colSums(apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))/N_k[k]
            
            xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
            psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
            rSinv_list <- mapply(S = S_list, 
                                 rik = as.list(r[k, ]), 
                                 FUN=function(S, rik){rik*solve(S)}, 
                                 SIMPLIFY=FALSE)
            rSinv_sum <- Reduce('+', rSinv_list)
            U_B[[k]] <- N_k[k]*d*solve(matrix(rowSums(mapply(x = xim, 
                                                             p = psim, 
                                                             rSinv = rSinv_list,
                                                             FUN=function(x,p,rSinv){
                                                                 v <- rbind(x, p)
                                                                 v%*%rSinv%*%t(v)  
                                                             }, SIMPLIFY=TRUE)), 
                                              nrow=2, byrow=FALSE))
            tryCatch(
                U_df[[k]] <- uniroot(function(nu0){(digamma_mv(x=nu0/2, p=d)
                                                    + 1/N_k[k]*sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                                                    - d*log(N_k[k]*nu0/2) 
                                                    + log(det(rSinv_sum))
                )}, lower = d+1, upper=1E9)$root, error=function(e){warning("Cluster too wide: inverse-Wishart degree of freedom very high")}
            )
            
            U_Sigma[[k]] <- N_k[k]*U_df[[k]]*solve(rSinv_sum)
        }
        
        loglik[i+1] <-sum(r*mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                                         U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                                         U_Sigma0 = U_Sigma, U_df0 = U_df))
        Q[i+1] <- (sum(r*kronecker(t(rep(1,ncol(r))), log(weights))) + loglik[i+1])
        
        
        
        cat("it ", i, ": Q = ", Q[i+1],"\n\n", sep="")
        if(abs(Q[i+1]-Q[i])<tol){break}
        
        plot(y=Q[2:(i+1)], x=c(1:i), 
             ylab="Q", xlab="Iteration", type="b", col="blue", pch=16)
        
    }
    
    plot(y=Q[2:(i+1)], x=c(1:i), 
         ylab="Q", xlab="it.", type="b", col="blue", pch=16)
    
    return(list("r"=r,
                "Q" = Q[2:(i+1)],
                "U_xi" = U_xi,
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df, 
                "U_Sigma" = U_Sigma,
                "weights"=weights))
    
}

#'MLE for sNiW distributed observations
#'
#'Maximum likelihood estimation of Normal inverse Wishart distributed observations
#'
#'@rdname MLE_skewT
#'
#'@export
#'
#'@examples
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(0.3, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 100
#'hyperG0$D_psi <- 100
#'hyperG0$nu <- 35
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'mle <- MLE_skewT(xi_list, psi_list, S_list)
#'mle

MLE_skewT <- function( xi_list, psi_list, S_list){
    
    
    N <- length(xi_list)
    d <- length(xi_list[[1]])
    
    if(length(psi_list) != N | length(S_list) != N){
        stop("Number of MCMC iterations not matching")
    }
    
    U_xi <- unlist(rowSums(sapply(xi_list, FUN="["))/N)
    U_psi <- unlist(rowSums(sapply(psi_list, FUN="["))/N)
    xim <- lapply(xi_list, function(x){x - U_xi})
    psim <- lapply(psi_list, function(x){x - U_psi})
    Sinv_list <- lapply(S_list, solve)
    Sinv_sum <- Reduce('+', Sinv_list)
    U_B <- N*d*solve(matrix(rowSums(mapply(x = xim, p = psim, Si = Sinv_list, FUN=function(x,p,Si){
        v <- rbind(x, p)
        tcrossprod(v%*%Si,v)  
    }, SIMPLIFY=TRUE)), 
    nrow=2, byrow=FALSE))
    tryCatch(
        U_df<- uniroot(function(nu0){(N/2*digamma_mv(x=nu0/2, p=d)
                                      + 1/2*sum(sapply(S_list, function(S){log(det(S))}))
                                      - N*d/2*log(N*nu0/2) 
                                      + N/2*log(det(Sinv_sum))
        )}, lower = d+1, upper=1E9)$root, error=function(e){warning("Cluster too wide: inverse-Wishart degree of freedom very high")}
    )
    
    U_Sigma <- N*U_df*solve(Sinv_sum)
    
    
    return(list("U_xi" = U_xi, 
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df, 
                "U_Sigma" = U_Sigma))
    
}
