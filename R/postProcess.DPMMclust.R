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
    w_list <- list()
    
    #m_final <- list()
    #S_final <- list()
    
    for(i in 1:length(x_invar$U_SS_list)){
        xi_list <- c(xi_list, sapply(x_invar$U_SS_list[[i]], "[", "xi"))
        psi_list <- c(psi_list, sapply(x_invar$U_SS_list[[i]], "[", "psi"))
        
        S_list <- c(S_list, sapply(x_invar$U_SS_list[[i]], "[", "S"))
        
        if(is.null(x_invar$U_SS_list[[1]][["weights"]])){
            #for compatibility with older DPMclust objects
            w_list <- c(w_list, x_invar$weights_list[[i]][unique(x_invar$mcmc_partitions[[i]])])
        }else{
            w_list <- c(w_list,sapply(x_invar$U_SS_list[[i]], "[", "weights"))
        }
        
        
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
#             param_post <- MAP_skewT_mmEM_weighted(xi_list, psi_list, S_list, 
#                                                   obsweight_list=w_list, 
#                                                   hyperG0 = x_invar$hyperG0, K=K, ...)
            param_post <- MAP_skewT_mmEM(xi_list, psi_list, S_list, 
                                                  hyperG0 = x_invar$hyperG0, K=K, ...)
        }
        else{
            param_post <- MLE_skewT(xi_list, psi_list, S_list, ...)
        }
    }
    else{
        stop("clust_distrib is not skewT\n other distrib nort implemented yet")
    }

    mle_g <- MLE_gamma(x_invar$alpha)

    parameters <- list()
    for (i in 1:length(param_post$U_xi)){
        parameters[[i]] <- list("b_xi" = param_post[["U_xi"]][[i]],
                                "b_psi" = param_post[["U_psi"]][[i]],
                                "B" = solve(param_post[["U_B"]][[i]]),
                                "lambda" = param_post[["U_Sigma"]][[i]],
                                "nu" = param_post[["U_df"]][[i]]
        )
    }
    
    return(list("parameters"=parameters, "weights"=param_post$weights,
                "alpha_param"=mle_g))
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
MLE_skewT_mmEM <- function( xi_list, psi_list, S_list, hyperG0, K, maxit=50, tol=1E-1, plot=TRUE){
    
    
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
    #Q <- numeric(maxit+1)
    #Q[1] <- -Inf
    
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
        #Q[i+1] <- (sum(r*kronecker(t(rep(1,ncol(r))), log(weights))) + loglik[i+1])
        
        
        cat("it ", i, ": loglik = ", loglik[i+1],"\n\n", sep="")
        if(abs(loglik[i+1]-loglik[i])<tol){break}
        if(plot){
            plot(y=loglik[2:(i+1)], x=c(1:i), 
                 ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
        }
    }
    if(plot){
        plot(y=loglik[2:(i+1)], x=c(1:i), 
             ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
    }
    
    return(list("r"=r,
                "loglik" = loglik[2:(i+1)],
                "U_xi" = U_xi,
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df, 
                "U_Sigma" = U_Sigma,
                "weights"=weights))
    
}

#'EM MAP for mixture of sNiW
#'
#'Maximum likelihood estimation of mixture of 
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@rdname MAP_skewT_mmEM
#'
#'@export
#'
#'@examples
#'set.seed(123)
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
#'w_list <- list()
#'
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'  w_list [[k]] <- 0.75
#'}
#'
#'
#'hyperG02 <- list()
#'hyperG02$b_xi <- c(-1, 2)
#'hyperG02$b_psi <- c(-0.1, 0.5)
#'hyperG02$kappa <- 0.001
#'hyperG02$D_xi <- 10
#'hyperG02$D_psi <- 10
#'hyperG02$nu <- 4
#'hyperG02$lambda <- 0.5*diag(2)
#'
#'for(k in 1001:2000){
#'  NNiW <- rNNiW(hyperG02, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'  w_list [[k]] <- 0.25
#'  
#'}
#'
#'library(lineprof)
#'map <- MAP_skewT_mmEM(xi_list, psi_list, S_list, hyperG0, K=4, tol=0.01)
#'map
#'
MAP_skewT_mmEM_vague <- function(xi_list, psi_list, S_list, hyperG0, K, maxit=50, tol=1E-1, plot=TRUE){
    
    
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
    
    
    #priors
    alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
    nu<- d+1
    lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=var))
    C <- diag(2)*1000
    L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=var)) 
          + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=var))
    )/2
    
    
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
    #Q <- numeric(maxit+1)
    #Q[1] <- -Inf
    
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
        weights  <- N_k/N #(N_k + alpha[k] - 1)/(N + sum(alpha) - K) #equivalent for alpha[k]=1
        
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
            U_B[[k]] <- (N_k[k]*d + 1)*solve(solve(C) + matrix(rowSums(mapply(x = xim, 
                                                                              p = psim, 
                                                                              rSinv = rSinv_list,
                                                                              FUN=function(x,p,rSinv){
                                                                                  v <- rbind(x, p)
                                                                                  v%*%rSinv%*%t(v)  
                                                                              }, SIMPLIFY=TRUE)), 
                                                               nrow=2, byrow=FALSE))
            U_df[[k]] <- try(uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                    + sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                                                    - N_k[k]*d*log(N_k[k]*nu0/2) 
                                                    + N_k[k]*log(det(rSinv_sum))
                                                    + 2
            )}, lower = d+1, upper=1E12)$root, TRUE)
            if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}
            
            
            
            U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
        }
        loglik[i+1] <- sum(log(apply(exp(mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                                                      U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                                                      U_Sigma0 = U_Sigma, U_df0 = U_df)), MARGIN=2, FUN=function(x){sum(x*weights)})))
        
        
        
        if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
            temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
            loglik[i+1] <- sum(temp_logliks)
        }
        
        cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
        cat("weights:", weights, "\n\n")
        if(abs(loglik[i+1]-loglik[i])<tol){break}
        
        if(plot){
            plot(y=loglik[2:(i+1)], x=c(1:i), 
                 ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
        }
    }
    
    if(plot){
        plot(y=loglik[2:(i+1)], x=c(1:i), 
             ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
    }
    
    return(list("r"=r,
                "loglik" = loglik[2:(i+1)],
                "U_xi" = U_xi,
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df,
                "U_Sigma" = U_Sigma,
                "weights"=weights))
    
}

#'@rdname MAP_skewT_mmEM
#'@export
MAP_skewT_mmEM<- function(xi_list, psi_list, S_list, hyperG0, K, maxit=50, tol=1E-1, plot=TRUE){
    
    
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
    
    
    #priors
    alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
    xi_p <- apply(sapply(xi_list, "["), MARGIN=1, FUN=mean)
    psi_p <- apply(sapply(psi_list, "["), MARGIN=1, FUN=mean)
    kappa0 <- 0.01
    nu<- d+1
    lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=var))
    C <- diag(2)*1000
    L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=var)) 
          + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=var))
    )/2
    
    
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
    #Q <- numeric(maxit+1)
    #Q[1] <- -Inf
    
    for(i in 1:maxit){
#         browser()
#         r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list, 
#                           U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
#                           U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))
        
        r <- mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                          U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                          U_Sigma0 = U_Sigma, U_df0 = U_df)
        r <- apply(X=r, MARGIN=2, FUN=function(x){x+log(weights)})
        r <- apply(X=r, MARGIN=2, FUN=function(x){x - log(sum(exp(x)))})
        r[which(is.infinite(r))] <- -Inf
        r <- exp(r)
        
        
        
        #M step
        N_k <- rowSums(r)
        weights  <- (N_k + alpha[k] - 1)/(N + sum(alpha) - K)
        
        for(k in 1:K){
            xi_m_k_xNk <- colSums(apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
            U_xi[[k]] <- (xi_m_k_xNk + kappa0/N*xi_p)/(N_k[k]  + kappa0/N)
            psi_m_k_xNk <- colSums(apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
            U_psi[[k]] <- (psi_m_k_xNk + kappa0/N*psi_p)/(N_k[k] + kappa0/N)
            
            xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
            psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
            xim0 <- U_xi[[k]] - xi_p
            psim0 <- U_xi[[k]] - psi_p
            
            Sinv_list <- lapply(S_list,solve)
            Sinv_sum <- Reduce('+', Sinv_list)
            rSinv_list <- mapply(Sinv = Sinv_list, 
                                 rik = as.list(r[k, ]), 
                                 FUN=function(Sinv, rik){rik*Sinv}, 
                                 SIMPLIFY=FALSE)
            rSinv_sum <- Reduce('+', rSinv_list)
            U_B[[k]] <- (N_k[k]*d + d + 1)*solve(solve(C) + matrix(rowSums(mapply(x = xim, 
                                                                                  p = psim, 
                                                                                  rSinv = rSinv_list,
                                                                                  FUN=function(x,p,rSinv){
                                                                                      v <- rbind(x, p)
                                                                                      v%*%rSinv%*%t(v)  
                                                                                  }, SIMPLIFY=TRUE)), 
                                                                   nrow=2, byrow=FALSE)
                                                 +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
            )
            U_df[[k]] <- try(uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                    + sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                                                    - N_k[k]*d*log(N_k[k]*nu0/2) 
                                                    + N_k[k]*log(det(rSinv_sum))
                                                    + 2
            )}, lower = d+1, upper=1E12)$root, TRUE)
            if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}
            
            
            
            U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
        }
        #        cat("df",unlist(U_df), "\n")
        #         loglik[i+1] <-sum(r*mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
        #                                          U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
        #                                          U_Sigma0 = U_Sigma, U_df0 = U_df))
        temp_logliks <- log(apply(exp(mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                                                   U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                                                   U_Sigma0 = U_Sigma, U_df0 = U_df)), MARGIN=2, FUN=function(x){sum(x*weights)}))
        loglik[i+1] <- sum(temp_logliks)
        
        #Q[i+1] <- (sum(r*kronecker(t(rep(1,ncol(r))), log(weights))) + loglik[i+1])
        
        
        
        if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
            temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
            loglik[i+1] <- sum(temp_logliks)
        }
        
        cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
        cat("weights:", weights, "\n\n")
        if(abs(loglik[i+1]-loglik[i])<tol){break}
        
        if(plot){
            plot(y=loglik[2:(i+1)], x=c(1:i), 
                 ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
        }
    }
    
    if(plot){
        plot(y=loglik[2:(i+1)], x=c(1:i), 
             ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
    }
    
    return(list("r"=r,
                "loglik" = loglik[2:(i+1)],
                "U_xi" = U_xi,
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df,
                "U_Sigma" = U_Sigma,
                "weights"=weights))
    
}

#'@rdname MAP_skewT_mmEM
#'@export
MAP_skewT_mmEM_weighted<- function(xi_list, psi_list, S_list, obsweight_list, hyperG0, K, maxit=50, tol=1E-1, plot=TRUE){
    
    
    pseudoN <- length(xi_list)
    N <- sum(rep(1,pseudoN)*unlist(obsweight_list))
    d <- length(hyperG0[[1]])
    
    if(length(psi_list) != pseudoN | length(S_list) != pseudoN | length(obsweight_list) != pseudoN){
        stop("Number of MCMC iterations not matching")
    }
    
    U_xi <- list() #matrix(0, nrow=d,ncol=K)
    U_psi <- list() #matrix(0, nrow=d,ncol=K)
    U_Sigma <- list() # array(dim=c(d,d,K))
    U_B <- list() #array(dim=c(2,2,K))
    U_df <- list() #numeric(K)
    
    #priors
    alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
    xi_p <- apply(sapply(xi_list, "["), MARGIN=1, FUN=mean)
    psi_p <- apply(sapply(psi_list, "["), MARGIN=1, FUN=mean)
    kappa0 <- 0.01
    nu<- d+1
    lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=var))
    C <- diag(2)*1000
    L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=var)) 
          + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=var))
    )/2
    
    
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
    #Q <- numeric(maxit+1)
    #Q[1] <- -Inf
    
    for(i in 1:maxit){
        
        r <- mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                          U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                          U_Sigma0 = U_Sigma, U_df0 = U_df)
        r <- apply(X=r, MARGIN=2, FUN=function(x){x + log(weights)})
        #r <- t(apply(X=r, MARGIN=1, FUN=function(x){x + log(unlist(obsweight_list))}))
        r <- apply(X=r, MARGIN=2, FUN=function(x){x - log(sum(exp(x)))})
        r[which(is.infinite(r))] <- -Inf
        r <- exp(r)
        r <- t(apply(X=r, MARGIN=1, FUN=function(x){x*unlist(obsweight_list)}))
        
        
        #M step
        N_k <- rowSums(r)
        weights  <- (N_k + alpha[k] - 1)/(N + sum(alpha) - K)
        
        for(k in 1:K){
            xi_m_k_xNk <- colSums(apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
            U_xi[[k]] <- (xi_m_k_xNk + kappa0/N*xi_p)/(N_k[k]  + kappa0/N)
            psi_m_k_xNk <- colSums(apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
            U_psi[[k]] <- (psi_m_k_xNk + kappa0/N*psi_p)/(N_k[k] + kappa0/N)
            
            xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
            psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
            xim0 <- U_xi[[k]] - xi_p
            psim0 <- U_xi[[k]] - psi_p
            
            Sinv_list <- lapply(S_list,solve)
            Sinv_sum <- Reduce('+', Sinv_list)
            rSinv_list <- mapply(Sinv = Sinv_list, 
                                 rik = as.list(r[k, ]), 
                                 FUN=function(Sinv, rik){rik*Sinv}, 
                                 SIMPLIFY=FALSE)
            rSinv_sum <- Reduce('+', rSinv_list)
            U_B[[k]] <- (N_k[k]*d + d + 1)*solve(solve(C) + matrix(rowSums(mapply(x = xim, 
                                                                                  p = psim, 
                                                                                  rSinv = rSinv_list,
                                                                                  FUN=function(x,p,rSinv){
                                                                                      v <- rbind(x, p)
                                                                                      v%*%rSinv%*%t(v)  
                                                                                  }, SIMPLIFY=TRUE)), 
                                                                   nrow=2, byrow=FALSE)
                                                 +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
            )
            U_df[[k]] <- try(uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                    + sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                                                    - N_k[k]*d*log(N_k[k]*nu0/2) 
                                                    + N_k[k]*log(det(rSinv_sum))
                                                    + 2
            )}, lower = d+1, upper=1E12)$root, TRUE)
            if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}
            
            
            
            U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
        }
        #        cat("df",unlist(U_df), "\n")
        #         loglik[i+1] <-sum(r*mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
        #                                          U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
        #                                          U_Sigma0 = U_Sigma, U_df0 = U_df))
        loglik[i+1] <- sum(log(apply(exp(mmsNiWlogpdf(U_xi = xi_list, U_psi = psi_list, U_Sigma = S_list, 
                                                      U_xi0 = U_xi, U_psi0 = U_psi, U_B0 =U_B,
                                                      U_Sigma0 = U_Sigma, U_df0 = U_df)), MARGIN=2, FUN=function(x){sum(x*weights)})))
        
        #Q[i+1] <- (sum(r*kronecker(t(rep(1,ncol(r))), log(weights))) + loglik[i+1])
        
        
        cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
        cat("weights:", weights, "\n\n")
        
        if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){browser()}
        if(abs(loglik[i+1]-loglik[i])<tol){break}
        
        if(plot){
            plot(y=loglik[2:(i+1)], x=c(1:i), 
                 ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
        }
    }
    
    if(plot){
        plot(y=loglik[2:(i+1)], x=c(1:i), 
             ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
    }
    
    return(list("r"=r,
                "loglik" = loglik[2:(i+1)],
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

MLE_skewT <- function( xi_list, psi_list, S_list, plot=TRUE){
    
    
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
    
    U_df<- try(uniroot(function(nu0){(N/2*digamma_mv(x=nu0/2, p=d)
                                      + 1/2*sum(sapply(S_list, function(S){log(det(S))}))
                                      - N*d/2*log(N*nu0/2) 
                                      + N/2*log(det(Sinv_sum))
    )}, lower = d+1, upper=1E9)$root, TRUE)
    if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}
    
    U_Sigma <- N*U_df*solve(Sinv_sum)
    
    
    return(list("U_xi" = U_xi, 
                "U_psi" = U_psi, 
                "U_B" = U_B, 
                "U_df" = U_df, 
                "U_Sigma" = U_Sigma))
    
}

#'MLE for Gamma distribution
#'
#'Maximum likelihood estimation of Gamma distributed observations 
#'distribution parameters
#'
#'
#'@export
#'
#'@examples
#'
#' g_list <- list()
#' for(i in 1:1000){
#'  g_list <- c(g_list, rgamma(1, shape=100, rate=5))
#' }
#' 
#' mle <- MLE_gamma(g_list)
#' mle
#'
MLE_gamma <- function(g){
    N <- length(g)
    
    a_mle <- try(uniroot(function(a){(N*mean(log(g))
                                     - N*digamma(a)
                                     - N*log(mean(g))
                                     + N*log(a)
    )}, lower = 0.000001, upper=1E9)$root, TRUE)
    if(inherits(a_mle, "try-error")){a_mle <- 0.0001;warning("unable to estimate a_mle properly")}    

    b_mle <- mean(g)/a_mle

    return(list("shape"=a_mle, "scale"=b_mle, "rate"=1/b_mle))
}
