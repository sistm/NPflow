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
#'partition of the \code{n} observations to compare to the point estimate.
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
#'@importFrom stats uniroot
#'
#'@export
#'
#'@importFrom gplots heatmap.2
#'
#'@seealso \code{\link{similarityMat}} \code{\link{summary.DPMMclust}}
#'
postProcess.DPMMclust <- function(x, burnin=0, thin=1, gs=NULL, lossFn="F-measure", K=10, ...){

  x_invar <- burn.DPMMclust(x, burnin = burnin, thin=thin)

  EM_init_nb_max <- 10
  elem <- which(lapply(x_invar$U_SS_list,FUN=length)==K)
  len <- length(elem)
  EM_init <- list()
  if(len>=EM_init_nb_max){
    EM_init_nb <- EM_init_nb_max
    randind <- elem[sample(1:len,EM_init_nb,replace=FALSE)]
    for (el in 1:EM_init_nb){
      EM_init[[el]] <- x_invar$U_SS_list[[randind[el]]]
    }
  }
  else{
    EM_init_nb <- len
    cpt <- 1
    for (el in elem){
      EM_init[[cpt]] <- x_invar$U_SS_list[[el]]
      cpt <- cpt+1
    }
  }


  if(x$clust_distrib=="skewt"){

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
    }

    mle_g <- MLE_gamma(x_invar$alpha)

    if(K>1){

      MAPprior <- x_invar$hyperG0
      #MAPprior$lambda <-10*MAPprior$lambda
      param_post_list <- list()

      chr_str <- paste(paste("/", as.character(EM_init_nb),sep=""), "computed",sep=" ")
      for (j in 1:EM_init_nb){
        param_post_list[[j]] <- MAP_skewT_mmEM(xi_list, psi_list, S_list,
                                               hyperG0 = MAPprior, K=K,
                                               init=EM_init[[j]],verbose=FALSE,...)

        cat("EM ", j,chr_str, "\n", sep="")
      }
      param_post <- param_post_list[[which.max(sapply(lapply(param_post_list, "[[", "loglik"), FUN=max))]]


      parameters <- list()
      for (i in 1:length(param_post$U_xi)){
        parameters[[i]] <- list("b_xi" = param_post[["U_xi"]][[i]],
                                "b_psi" = param_post[["U_psi"]][[i]],
                                "B" = param_post[["U_B"]][[i]],
                                "lambda" = param_post[["U_Sigma"]][[i]],
                                "nu" = param_post[["U_df"]][[i]]
        )
      }
    }
    else{
      param_post <- MLE_skewT(xi_list, psi_list, S_list, ...)
      parameters <- list()
      parameters[[1]] <- list("b_xi" = param_post[["U_xi"]],
                              "b_psi" = param_post[["U_psi"]],
                              "B" = param_post[["U_B"]],
                              "lambda" = param_post[["U_Sigma"]],
                              "nu" = param_post[["U_df"]]
      )
      param_post$weights <- 1
    }


  }else if (x$clust_distrib=="gaussian"){

    mle_g <- MLE_gamma(x_invar$alpha)

    mu_list <- list()
    S_list <- list()
    w_list <- list()


    for(i in 1:length(x_invar$U_SS_list)){
      mu_list <- c(mu_list, sapply(x_invar$U_SS_list[[i]], "[", "mu"))
      S_list <- c(S_list, sapply(x_invar$U_SS_list[[i]], "[", "S"))

      if(is.null(x_invar$U_SS_list[[1]][["weights"]])){
        #for compatibility with older DPMclust objects
        w_list <- c(w_list, x_invar$weights_list[[i]][unique(x_invar$mcmc_partitions[[i]])])
      }else{
        w_list <- c(w_list,sapply(x_invar$U_SS_list[[i]], "[", "weights"))
      }
    }

    param_post_list <- list()
    for (j in 1:EM_init_nb){

      param_post_list[[j]] <- MLE_N_mmEM(mu_list, S_list, x_invar$hyperG0, K, maxit=100, tol=1E-1, plot=TRUE)
      cat("EM ", j, "/10 computed", "\n", sep="")
    }
    param_post <- param_post_list[[which.max(sapply(lapply(param_post_list, "[[", "loglik"), FUN=max))]]


    parameters <- list()
    for (i in 1:length(param_post$U_mu)){
      parameters[[i]] <- list("mu" = as.vector(param_post[["U_mu"]][[i]]),
                              "kappa" = param_post[["U_kappa"]][[i]],
                              "lambda" = param_post[["U_lambda"]][[i]],
                              "nu" = param_post[["U_nu"]][[i]]
      )
    }

  }
  else {stop("clust_distrib is neither 'skewt' nor 'gaussian'\n other distributions are not implemented yet")}

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
#'@importFrom stats uniroot
#'
#'@importFrom graphics plot
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
#'hyperG0$nu <- 3
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
MLE_skewT_mmEM <- function( xi_list, psi_list, S_list, hyperG0, K, init=NULL,maxit=100, tol=1E-1, plot=TRUE,verbose=TRUE){


  N <- length(xi_list)
  d <- length(xi_list[[1]])

  if(length(psi_list) != N | length(S_list) != N){
    stop("Number of MCMC iterations not matching")
  }

  U_xi <- list() #matrix(0, nrow=d,ncol=K)
  U_psi <- list() #matrix(0, nrow=d,ncol=K)
  U_Sigma <- list() # array(dim=c(d,d,K))
  U_B <- list() #array(dim=c(2,2,K))
  U_df <- list() #numeric(K)


  #initialisation
  if(is.null(init)){
    for(k in 1:K){
      #sampling the cluster parameters
      NNiW <- rNNiW(hyperG0, diagVar=FALSE)
      U_xi[[k]] <- NNiW[["xi"]]
      U_psi[[k]] <- NNiW[["psi"]]
      U_Sigma[[k]] <- NNiW[["S"]]
      U_B[[k]] <- diag(100, 2)
      U_df[[k]] <- d+1
    }
  }
  else{
    for(k in 1:K){
      #cluster parameters
      U_xi[[k]] <- init[[k]]$b_xi
      U_psi[[k]] <- init[[k]]$b_psi
      U_Sigma[[k]] <-  init[[k]]$lambda
      U_B[[k]] <- init[[k]]$B
      U_df[[k]] <- init[[k]]$nu
    }
  }

  weights <- rep(1/K, K)
  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))



  #Q <- numeric(maxit+1)
  #Q[1] <- -Inf


  r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                  U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                  U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

  logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
  r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
  rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
  r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))
  N_k<-rowSums(r)




  for(i in 1:maxit){

    #E step
    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))



    #M step
    N_k <- rowSums(r)
    weights  <- N_k/N
    #     cat("weights:", weights, "\n")


    for(k in 1:K){

      rSinv_list <- mapply(S = S_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*solve(S)},
                           SIMPLIFY=FALSE)

      rSinv_sum <- Reduce('+', rSinv_list)


      r_xi<-apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_psi<-apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})

      r_xi_psi<-cbind(r_xi,r_psi)


      U_xi_U_psi<-Reduce('+',lapply(1:N,
                                    function(i,m) {t(matrix(r_xi_psi[i,],d,d))%*%m[[i]]},
                                    m=rSinv_list))%*% solve(rSinv_sum)


      U_xi[[k]] <- U_xi_U_psi[1,]
      U_psi[[k]]<- U_xi_U_psi[2,]
      xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
      psim <- lapply(psi_list, function(x){x - U_psi[[k]]})

      #       U_B[[k]] <- 1/(N_k[k]*d)*(matrix(rowSums(mapply(x = xim,
      #                                                       p = psim,
      #                                                       rSinv = rSinv_list,
      #                                                       FUN=function(x,p,rSinv){
      #                                                         v <- rbind(x, p)
      #                                                         v%*%rSinv%*%t(v)
      #                                                       }, SIMPLIFY=TRUE)),
      #                                        nrow=2, byrow=FALSE))

      U_B[[k]] <- 1/(N_k[k]*d)*(diag(diag((matrix(rowSums(mapply(x = xim,
                                                                 p = psim,
                                                                 rSinv = rSinv_list,
                                                                 FUN=function(x,p,rSinv){
                                                                   v <- rbind(x, p)
                                                                   v%*%rSinv%*%t(v)
                                                                 }, SIMPLIFY=TRUE)),
                                                  nrow=2, byrow=FALSE)))))




      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)

      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                              - N_k[k]*d*log(N_k[k]*nu0/2)
                                              + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}

      U_Sigma[[k]] <- N_k[k]*U_df[[k]]*solve(rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


    #     cat("it ", i, ": loglik = ", loglik[i+1],"\n\n", sep="")
    if(abs(loglik[i+1]-loglik[i])<tol){break}
    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }

    if(plot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
           ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }
  if(plot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
         ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
  }

  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}


#'@rdname MAP_skewT_mmEM
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
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
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(1, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.1
#'hyperG0$D_xi <- 1
#'hyperG0$D_psi <- 1
#'hyperG0$nu <- 2
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
#'hyperG02$kappa <- 0.1
#'hyperG02$D_xi <- 1
#'hyperG02$D_psi <- 1
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
#'mle <- MLE_skewT_mmEM(xi_list, psi_list, S_list, hyperG0, K=2, tol=0.1)
#'mle
#'
MAP_skewT_mmEM_vague <- function(xi_list, psi_list, S_list,
                                 hyperG0, K, maxit=100, tol=1E-1, plot=TRUE){


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
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters
    NNiW <- rNNiW(hyperG0, diagVar=FALSE)
    U_xi[[k]] <- NNiW[["xi"]]
    U_psi[[k]] <- NNiW[["psi"]]
    U_Sigma[[k]] <- NNiW[["S"]]
    U_B[[k]] <- diag(100, 2)
    U_df[[k]] <- d+1
  }


  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))
  #Q <- numeric(maxit+1)
  #Q[1] <- -Inf

  for(i in 1:maxit){
    #browser()

    #E step
    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))



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
      U_B[[k]] <- 1/(N_k[k]*d + 1)*(solve(C) + matrix(rowSums(mapply(x = xim,
                                                                     p = psim,
                                                                     rSinv = rSinv_list,
                                                                     FUN=function(x,p,rSinv){
                                                                       v <- rbind(x, p)
                                                                       v%*%rSinv%*%t(v)
                                                                     }, SIMPLIFY=TRUE)),
                                                      nrow=2, byrow=FALSE))
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)
      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                              - N_k[k]*d*log(N_k[k]*nu0/2)
                                              + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}


      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))



    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
      loglik[i+1] <- sum(temp_logliks)
    }

    cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
    cat("weights:", weights, "\n\n")
    if(abs(loglik[i+1]-loglik[i])<tol){break}

    if(plot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
           ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }

  if(plot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
         ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
  }

  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}

#'EM MAP for mixture of sNiW
#'
#'Maximum A Posteriori (MAP) estimation of mixture of
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@param xi_list a list of length \code{n}, each element is a vector of size \code{d}
#'containing the argument \code{xi} of the corresponding allocated cluster.
#'
#'@param psi_list a list of length \code{n}, each element is a vector of size \code{d}
#'containing the argument \code{psi} of the corresponding allocated cluster.
#'
#'@param S_list a list of length \code{n}, each element is a matrix of size \code{d x d}
#'containing the argument \code{S} of the corresponding allocated cluster.
#'
#'@param hyperG0 prior mixing distribution.
#'
#'@param K integer giving the number of mixture components.
#'
#'@param maxit integer giving the maximum number of iteration for the EM algorithm. Default is 100.
#'
#'@param tol real number giving the tolerance for the stopping of the EM algorithm.
#'Default is \code{0.1}.
#'
#'@param doPlot a logical flag indicating wether the algorithm progresison should be plotted.
#'Default is \code{TRUE}.
#'
#'@rdname MAP_skewT_mmEM
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
#'
#'@export
MAP_skewT_mmEM<- function(xi_list, psi_list, S_list, hyperG0, init=NULL, K, maxit=100, tol=1E-1, plot=TRUE, verbose=TRUE){



  N <- length(xi_list)
  d <- length(xi_list[[1]])

  if(length(psi_list) != N | length(S_list) != N){
    stop("Number of MCMC iterations not matching")
  }

  U_xi <- list() #matrix(0, nrow=d,ncol=K)
  U_psi <- list() #matrix(0, nrow=d,ncol=K)
  U_Sigma <- list() # array(dim=c(d,d,K))
  U_B <- list() #array(dim=c(2,2,K))
  U_df <- list() #numeric(K)


  #priors
  #   alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  xi_p <- apply(sapply(xi_list, "["), MARGIN=1, FUN=mean)
  psi_p <- apply(sapply(psi_list, "["), MARGIN=1, FUN=mean)
  kappa0 <- 0.01
  nu<- d+1
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation


  if(is.null(init)){
    for(k in 1:K){
      #sampling the cluster parameters
      NNiW <- rNNiW(hyperG0, diagVar=FALSE)
      U_xi[[k]] <- NNiW[["xi"]]
      U_psi[[k]] <- NNiW[["psi"]]
      U_Sigma[[k]] <- NNiW[["S"]]
      U_B[[k]] <- diag(100, 2)
      U_df[[k]] <- d+1
    }
  }
  else{
    for(k in 1:K){
      #cluster parameters
      U_xi[[k]] <- init[[k]]$b_xi
      U_psi[[k]] <- init[[k]]$b_psi
      U_Sigma[[k]] <- init[[k]]$lambda
      U_B[[k]] <- init[[k]]$B
      U_df[[k]] <- init[[k]]$nu
    }
  }


  K<-K
  #priors
  alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  weights <- rep(1/K, K)

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


  for(i in 1:maxit){

    #E step

    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))


    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


    #M step
    N_k <- rowSums(r)
    weights  <- (N_k + alpha[k] - 1)/(N + sum(alpha) - K)

    for(k in 1:K){


      Sinv_list <- mapply(S = S_list,
                          FUN=function(S){solve(S)},
                          SIMPLIFY=FALSE)
      Sinv_sum <- Reduce('+', Sinv_list)
      rSinv_list <- mapply(S = Sinv_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*S},
                           SIMPLIFY=FALSE)
      rSinv_sum <- Reduce('+', rSinv_list)
      r_xi<-apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_psi<-apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_xi_psi<-cbind(r_xi,r_psi)
      U_xi_U_psi<-((kappa0/N*t(matrix(cbind(xi_p,psi_p),d,d))%*%Sinv_sum)+Reduce('+',lapply(1:N,
                                                                                            function(i,m) {t(matrix(r_xi_psi[i,],d,d))%*%m[[i]]},
                                                                                            m=rSinv_list)))%*% solve(kappa0/N*Sinv_sum+rSinv_sum)

      U_xi[[k]] <- U_xi_U_psi[1,]
      U_psi[[k]]<- U_xi_U_psi[2,]

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
      U_B[[k]] <- 1/(N_k[k]*d + d + 1)*(solve(C) + matrix(rowSums(mapply(x = xim,  #1/
                                                                         p = psim,
                                                                         rSinv = rSinv_list,
                                                                         FUN=function(x,p,rSinv){
                                                                           v <- rbind(x, p)
                                                                           v%*%rSinv%*%t(v)
                                                                         }, SIMPLIFY=TRUE)),
                                                          nrow=2, byrow=FALSE)
                                        +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
      )
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)

      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                              - N_k[k]*d*log(N_k[k]*nu0/2)
                                              + const_nu0_uniroot
      )}, lower = d, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}



      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }
    #        cat("df",unlist(U_df), "\n")

    temp_logliks <- log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                     U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                     U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE)))
    loglik[i+1] <- sum(temp_logliks)


    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      #       browser()
      temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
      loglik[i+1] <- sum(temp_logliks)
    }

    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }

    if(plot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
           ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }

    if(abs(loglik[i+1]-loglik[i])<tol){break}
  }
  U_B <- lapply(U_B,function(i)(diag(diag(i))))
  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}


#'@rdname MAP_skewT_mmEM
#'
#'@param obsweight_list a list of length \code{n} where each element is a vector of weights for
#'each sampled cluster at each MCMC iterations.
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
#'
#'@export
MAP_skewT_mmEM_weighted<- function(xi_list, psi_list, S_list, obsweight_list,
                                   hyperG0, K, maxit=100, tol=1E-1, doPlot=TRUE){


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
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters
    NNiW <- rNNiW(hyperG0, diagVar=FALSE)
    U_xi[[k]] <- NNiW[["xi"]]
    U_psi[[k]] <- NNiW[["psi"]]
    U_Sigma[[k]] <- NNiW[["S"]]
    U_B[[k]] <- diag(100, 2)
    U_df[[k]] <- d+1
  }

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))

  for(i in 1:maxit){

    #E step
    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


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
      U_B[[k]] <- 1/(N_k[k]*d + d + 1)*solve(solve(C) + matrix(rowSums(mapply(x = xim,
                                                                              p = psim,
                                                                              rSinv = rSinv_list,
                                                                              FUN=function(x,p,rSinv){
                                                                                v <- rbind(x, p)
                                                                                v%*%rSinv%*%t(v)
                                                                              }, SIMPLIFY=TRUE)),
                                                               nrow=2, byrow=FALSE)
                                             +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
      )
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)
      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                              - N_k[k]*d*log(N_k[k]*nu0/2)
                                              + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}


      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


    cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
    cat("weights:", weights, "\n\n")

    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      #browser()
    }
    if(abs(loglik[i+1]-loglik[i])<tol){break}

    if(plot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
           ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }

  if(plot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
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
#'@importFrom stats uniroot
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

  Sinv_list <- lapply(S_list, solve)
  Sinv_sum <- Reduce('+', Sinv_list)

  U_xi_U_psi<-Reduce('+',lapply(1:N,
                                function(i,m) {rbind(xi_list[[i]],psi_list[[i]])%*%m[[i]]},
                                m=Sinv_list))%*% solve(Sinv_sum)
  U_xi<-U_xi_U_psi[1,]
  U_psi<-U_xi_U_psi[2,]

  xim <- lapply(xi_list, function(x){x - U_xi})
  psim <- lapply(psi_list, function(x){x - U_psi})

  U_B <- N*d*solve(matrix(rowSums(mapply(x = xim, p = psim, Si = Sinv_list, FUN=function(x,p,Si){
    v <- rbind(x, p)
    tcrossprod(v%*%Si,v)
  }, SIMPLIFY=TRUE)),
  nrow=2, byrow=FALSE))

  U_df<- try(stats::uniroot(function(nu0){(N/2*digamma_mv(x=nu0/2, p=d)
                                    + 1/2*sum(sapply(S_list, function(S){log(det(S))}))
                                    - N*d/2*log(N*nu0/2)
                                    + N/2*log(det(Sinv_sum))
  )}, lower = d+1, upper=1E9)$root, TRUE)
  if(inherits(U_df, "try-error")){U_df <- d+1}

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
#'@param g a list of Gamma distributed observation.
#'
#'@importFrom stats uniroot
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

  a_mle <- try(stats::uniroot(function(a){(N*mean(log(g))
                                    - N*digamma(a)
                                    - N*log(mean(g))
                                    + N*log(a)
  )}, lower = 0.000001, upper=1E9)$root, TRUE)
  if(inherits(a_mle, "try-error")){
    #either
    a_mle <- 0.0001
    b_mle <- 0.0001
    warning("Unable to estimate Gamma hyperpriors properly
            (this can happen when only a few clusters are fitted).
            => Non informative values are returned instead")
  }else{
    b_mle <- mean(g)/a_mle
  }

  return(list("shape"=a_mle, "scale"=b_mle, "rate"=1/b_mle))
}




#'Maximum likelihood estimation of mixture of
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@rdname MLE_N_mmEM
#'
#'@importFrom stats uniroot
#'
#'@importFrom graphics plot
#'
#'@export
#'
#'@examples
#'set.seed(123)
#' U_mu <- list()
#' U_Sigma <- list()
#' U_nu<-list()
#' U_kappa<-list()
#'
#' d <- 2
#' hyperG0 <- list()
#' hyperG0[["mu"]] <- rep(1,d)
#' hyperG0[["kappa"]] <- 0.01
#' hyperG0[["nu"]] <- d+1
#' hyperG0[["lambda"]] <- diag(d)
#'
#' for(k in 1:1000){
#'
#'   NiW <- rNiW(hyperG0, diagVar=FALSE)
#'   U_mu[[k]] <-NiW[["mu"]]
#'   U_Sigma[[k]] <-NiW[["S"]]
#' }
#'
#'
#' hyperG02 <- list()
#' hyperG02[["mu"]] <- rep(2,d)
#' hyperG02[["kappa"]] <- 1
#' hyperG02[["nu"]] <- d+10
#' hyperG02[["lambda"]] <- diag(d)/10
#'
#' for(k in 1001:2000){
#'
#'   NiW <- rNiW(hyperG02, diagVar=FALSE)
#'   U_mu[[k]] <-NiW[["mu"]]
#'   U_Sigma[[k]] <-NiW[["S"]]
#' }
#'
#'
#' mle <- MLE_N_mmEM( U_mu, U_Sigma, hyperG0, K=2)
#'
#' hyperG0[["mu"]]
#' hyperG02[["mu"]]
#' mle$U_mu
#'
#' hyperG0[["lambda"]]
#' hyperG02[["lambda"]]
#' mle$U_lambda
#'
#' hyperG0[["nu"]]
#' hyperG02[["nu"]]
#' mle$U_nu
#'
#' hyperG0[["kappa"]]
#' hyperG02[["kappa"]]
#' mle$U_kappa


MLE_N_mmEM <- function( mu_list, S_list, hyperG0, K, maxit=100, tol=1e-1, plot=TRUE){

  N <- length(mu_list)
  d <- length(hyperG0[[1]])

  if(length(mu_list) != N | length(S_list) != N){
    stop("Number of MCMC iterations not matching")
  }

  U_mu <- list()
  U_Sigma <- list()
  U_nu<-list()
  U_kappa<-list()

  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters

    NiW <- rNiW(hyperG0, diagVar=FALSE)
    U_mu[[k]] <-NiW[["mu"]]
    U_Sigma[[k]] <-NiW[["S"]]
    U_nu[[k]]<-d+1
    U_kappa[[k]]<-0.01
  }

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
                                                  U_mu0=sapply(U_mu, "["),
                                                  U_kappa0=sapply(U_kappa, "["),
                                                  U_nu0=sapply(U_nu, "["),
                                                  U_lambda0=U_Sigma,Log=FALSE))))
  loglik_w <- loglik

  for(i in 1:maxit){

    #E step
    #     r <- mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
    #                      U_mu0=sapply(U_mu, "["),
    #                      U_kappa0=sapply(U_kappa, "["),
    #                      U_nu0=sapply(U_nu, "["),
    #                      U_lambda0=U_Sigma,Log=TRUE)


    r <- mmNiWpdfC(Mu=sapply(mu_list, "["), Sigma=S_list,
                   U_Mu0=sapply(U_mu, "["),
                   U_Kappa0=sapply(U_kappa, "["),
                   U_Nu0=sapply(U_nu, "["),
                   U_Sigma0=U_Sigma,Log=TRUE)


    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


    #M step
    N_k <- rowSums(r)
    if(any(is.nan(N_k)) | any(round(N_k)==1) | any(!N_k) |any(rowSums((r==0)*1)==length(r[1,])-1)){
      # the above line is because sometimes for a given group, there is only one (or zero) loglikelihood superior to the likelihoods
      # computed for the other groups
      for(k in 1:K){
        #sampling the cluster parameters

        NiW <- rNiW(hyperG0, diagVar=FALSE)
        U_mu[[k]] <-NiW[["mu"]]
        U_Sigma[[k]] <-NiW[["S"]]
        U_nu[[k]]<-d+1
        U_kappa[[k]]<-0.01
      }
      maxit<-maxit+1
      next
    }
    weights  <- N_k/N
    #cat("weights:", weights, "\n")

    for(k in 1:K){

      rSinv_list <- mapply(S = S_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*solve(S)},
                           SIMPLIFY=FALSE)

      rSinv_sum <- Reduce('+', rSinv_list)


      U_mu[[k]]<-Reduce('+',lapply(1:N,
                                   function(i,m) {r[k,i]*mu_list[[i]]%*%m[[i]]},
                                   m=rSinv_list))%*% solve(rSinv_sum)


      mu_m <- lapply(mu_list, function(x){x - as.vector(U_mu[[k]])})

      U_kappa[[k]] <- (N_k[k]*d)/sum(mapply(x = mu_m, rSinv = rSinv_list,
                                            FUN=function(x,rSinv){
                                              t(x)%*%rSinv%*%x
                                            }, SIMPLIFY=TRUE))

      const_nu0_uniroot <-(sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                           + N_k[k]*log(det(rSinv_sum)))


      max_upper<-1e12
      U_nu[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                              - N_k[k]*d*log(N_k[k]*nu0/2)
                                              + const_nu0_uniroot
      )}, lower = d+1, upper=max_upper)$root, TRUE)


      if(inherits(U_nu[[k]], "try-error")){U_nu[[k]] <- d+1}
      U_nu[[k]]<-round(U_nu[[k]])  # warning issued in wishrnd: n isn't an integer anymore


      U_Sigma[[k]] <- N_k[k]*U_nu[[k]]*solve(rSinv_sum)
    }



    loglik[i+1] <- sum(log(colSums(weights%*%mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
                                                      U_mu0=sapply(U_mu, "["),
                                                      U_kappa0=sapply(U_kappa, "["),
                                                      U_nu0=sapply(U_nu, "["),
                                                      U_lambda0=U_Sigma,Log=FALSE))))


    if (is.nan(loglik[i+1])){break}
    #cat("it ", i, ": loglik = ", loglik[i+1],"\n\n", sep="")
    if(abs(loglik[i+1]-loglik[i])<tol){break}
    #if (loglik[i+1]<loglik[i]){browser()}
    if(plot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
           ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }
  if(plot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
         ylab="Log-likelihood", xlab="it.", type="b", col="blue", pch=16)
  }

  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_mu" = U_mu,
              "U_kappa"=U_kappa,
              "U_nu" = U_nu,
              "U_lambda" = U_Sigma,
              "weights"=weights))

}
