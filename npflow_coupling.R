obsdim <- 2
nobservations <- 50
z <- matrix(rnorm(nobservations*obsdim), ncol = nobservations, nrow = obsdim)
diagVar=TRUE
use_variance_hyperprior <- TRUE
a=0.0001
b=0.0001
hyperG0 <- list()
hyperG0[["mu"]] <- rep(0,obsdim)
hyperG0[["kappa"]] <- 0.001
hyperG0[["nu"]] <- obsdim+2
hyperG0[["lambda"]] <- diag(obsdim)/10
nbclust_init <- 5
p <- nrow(z)
n <- ncol(z)

## init Gibbs
rinit <- function(){
  U_mu <- matrix(0, nrow=p, ncol=n)
  U_Sigma = array(0, dim=c(p, p, n))
  U_SS <- list()
  m <- numeric(n) # number of obs in each clusters
  c <- numeric(n)
  c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
  for (k in unique(c)){
    obs_k <- which(c==k)
    #cat("cluster ", k, ":\n")
    U_SS[[k]] <- NPflow:::update_SS(z=z[, obs_k,drop=FALSE], S=hyperG0)
    NiW <- NPflow:::rNiW(U_SS[[k]], diagVar=diagVar)
    
    U_mu[, k] <- NiW[["mu"]]
    U_SS[[k]][["mu"]] <- NiW[["mu"]]
    
    U_Sigma[, , k] <- NiW[["S"]]
    U_SS[[k]][["S"]] <- NiW[["S"]]
    
    m[k] <- length(obs_k)
    # U_SS[[k]][["weight"]] <- m[k]/n
  }
  
  alpha <- log(n)
  ## state of the markov chain:
  ## U_mu, U_Sigma, U_SS, c, alpha
  return(list(alpha = alpha, c = c, U_mu = U_mu, U_Sigma = U_Sigma, U_SS = U_SS, m = m))
}  

## Markov kernel
## state -> state 

single_kernel <- function(state){
  alpha <- state$alpha
  c <- state$c
  U_mu <- state$U_mu
  U_Sigma <- state$U_Sigma
  U_SS <- state$U_SS
  nbClust <- length(unique(c))
  # update alpha
  alpha <- NPflow:::sample_alpha(alpha_old=alpha, n=n, K=nbClust, a=a, b=b)
  # slice sampling update (c)
  slice <- NPflow:::sliceSampler_N(c=c, m=m, alpha=alpha, z=z, hyperG0=hyperG0,
                                   U_mu=U_mu, U_Sigma=U_Sigma, diagVar=diagVar)
  m <- slice[["m"]]
  c <- slice[["c"]]
  # weights_list[[i]] <- slice[["weights"]]
  U_mu<-slice[["U_mu"]]
  U_Sigma<-slice[["U_Sigma"]]
  # Update cluster locations
  fullCl <- which(m!=0)
  for(j in fullCl){
    obs_j <- which(c==j)
    #cat("cluster ", j, ":\n")
    if (use_variance_hyperprior){
      U_SS[[j]] <- NPflow:::update_SS(z=z[, obs_j, drop=FALSE], S=hyperG0, hyperprior = list("Sigma"=U_Sigma[,,j]))
    } else{
      U_SS[[j]] <- NPflow:::update_SS(z=z[, obs_j, drop=FALSE], S=hyperG0)
    }
    NiW <- NPflow:::rNiW(U_SS[[j]], diagVar=diagVar)
    U_mu[, j] <- NiW[["mu"]]
    U_SS[[j]][["mu"]] <- NiW[["mu"]]
    
    U_Sigma[, , j] <- NiW[["S"]]
    U_SS[[j]][["S"]] <- NiW[["S"]]
    
    # U_SS[[j]][["weight"]] <- weights_list[[i]][j]
    #cat("sampled S =", NiW[["S"]], "\n\n\n")
  }
  return(list(alpha = alpha, c = c, U_mu = U_mu, U_Sigma = U_Sigma, U_SS = U_SS, m = m))
}

# state <- rinit()
# nmcmc <- 1e2
# for (imcmc in 1:nmcmc){
#   state <- single_kernel(state)
#   # + store 
# }

## Coupled Markov kernel
## state1, state2 -> state1, state2 
get_max_coupling <- function(rp, dp, rq, dq){
  x <- rp()
  if (dp(x) + log(runif(1)) < dq(x)){
    return(list(xy = c(x,x), identical = TRUE))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rq()
      reject <- (dq(y) + log(runif(1)) < dp(y))
    }
    return(list(xy = c(x,y), identical = FALSE))
  }
}

gamma_max_coupling <- function(shape1, shape2, scale1, scale2){
  get_max_coupling(rp = function() rgamma(1, shape1, scale = scale1),
                   dp = function(x) dgamma(x, shape1, scale = scale1, log = TRUE),
                   rq = function() rgamma(1, shape2, scale = scale2),
                   dq = function(x) dgamma(x, shape2, scale = scale2, log = TRUE))  
}

sample_coupled_alpha <- function(alpha_old1, alpha_old2, n, K1, K2, a=0.0001, b=0.0001){
  
  # if(b > 0){
  # Sample scale factor in Dirichlet Process
  # x1 <- stats::rbeta(n=1, shape1=alpha_old1 + 1, shape2=n)
  # x2 <- stats::rbeta(n=1, shape1=alpha_old2 + 1, shape2=n)
  rp <- function() rbeta(n=1, shape1=alpha_old1 + 1, shape2=n)
  dp <- function(x) dbeta(x=x, shape1=alpha_old1 + 1, shape2=n, log = TRUE)
  rq <- function() rbeta(n=1, shape1=alpha_old2 + 1, shape2=n)
  dq <- function(x) dbeta(x=x, shape1=alpha_old2 + 1, shape2=n, log = TRUE)
  X1X2 <- get_max_coupling(rp, dp, rq, dq)
  x1 <- X1X2$xy[1]
  x2 <- X1X2$xy[2]
  temp1 <- (a+K1-1) / (n*(b-log(x1)))
  pi1 <- temp1/(1+temp1)
  temp2 <- (a+K2-1) / (n*(b-log(x2)))
  pi2 <- temp2/(1+temp2)
  u <- stats::runif(1)
  if (u < pi1 && u < pi2){
    # alpha_new1  <-  stats::rgamma(n=1, shape=a + K1, scale=1/(b - log(x1)))
    # alpha_new2  <-  stats::rgamma(n=1, shape=a + K2, scale=1/(b - log(x2)))
    alpha1alpha2 <- gamma_max_coupling(shape1 = a + K1, shape2 = a + K2, scale1 = 1/(b - log(x1)), scale2 = 1/(b - log(x2)))
    alpha_new1 <- alpha1alpha2$xy[1]
    alpha_new2 <- alpha1alpha2$xy[2]
  }
  if (u < pi1 && u > pi2){
    alpha1alpha2 <- gamma_max_coupling(shape1 = a + K1, shape2 = a + K2 - 1, scale1 = 1/(b - log(x1)), scale2 = 1/(b - log(x2)))
    alpha_new1 <- alpha1alpha2$xy[1]
    alpha_new2 <- alpha1alpha2$xy[2]
  }
  if (u > pi1 && u < pi2){
    alpha1alpha2 <- gamma_max_coupling(shape1 = a + K1 - 1, shape2 = a + K2, scale1 = 1/(b - log(x1)), scale2 = 1/(b - log(x2)))
    alpha_new1 <- alpha1alpha2$xy[1]
    alpha_new2 <- alpha1alpha2$xy[2]
    
  }
  if (u > pi1 && u > pi2){
    alpha1alpha2 <- gamma_max_coupling(shape1 = a + K1 - 1, shape2 = a + K2 -1, scale1 = 1/(b - log(x1)), scale2 = 1/(b - log(x2)))
    alpha_new1 <- alpha1alpha2$xy[1]
    alpha_new2 <- alpha1alpha2$xy[2]
  }
  return(c(alpha_new1, alpha_new2))
  # }else if(b==0){
  #   alpha_new=a
  # }else{
  #   print("b cannot be negative")
  # }
}

state1 <- rinit()
state1 <- single_kernel(state1)
state2 <- rinit()
state2 <- single_kernel(state2)

m1 = state1$m
m2 = state2$m
c1 = state1$c
c2 = state2$c
alpha1 = state1$alpha
alpha2 = state2$alpha
U_mu1 = state1$U_mu
U_mu2 = state2$U_mu
U_Sigma1 = state1$U_Sigma
U_Sigma2 = state2$U_Sigma

coupled_sliceSampler_N <- function(c1, c2, m1, m2, alpha1, alpha2, z, hyperG0, U_mu1, U_mu2, U_Sigma1, U_Sigma2, diagVar){
  
  maxCl1 <- length(m1) # maximum number of clusters
  ind1 <- which(m1!=0) # indexes of non empty clusters
  maxCl2 <- length(m2) # maximum number of clusters
  ind2 <- which(m2!=0) # indexes of non empty clusters
  
  # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
  # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1)    # and sample the rest of the weigth for potential new clusters:
  # temp_{K+1} ~ Gamma(alpha, 1)
  # then renormalise temp
  w1 <- numeric(maxCl1)
  w2 <- numeric(maxCl2)
  shapes1 <- c(m1[ind1], alpha1)
  shapes2 <- c(m2[ind2], alpha2)
  temp1 <- numeric(length(shapes1))
  temp2 <- numeric(length(shapes2))
  for (i in 1:min(length(shapes1), length(shapes2))){
    xy_ <- gamma_max_coupling(shape1 = shapes1[i], shape2 = shapes2[i], scale1 = 1, scale2 = 1)
    temp1[i] <- xy_$xy[1]
    temp2[i] <- xy_$xy[2]
  }
  if (length(shapes1)!= length(shapes2)){
    if (length(shapes1) < length(shapes2)){
      temp2[length(shapes1):length(shapes2)] <- rgamma(n = length(shapes2)-length(shapes1), shape = shapes1[length(shapes1):length(shapes2)], scale = 1)
    } else {
      temp1[length(shapes2):length(shapes1)] <- rgamma(n = length(shapes1)-length(shapes2), shape = shapes2[length(shapes2):length(shapes1)], scale = 1)
    }
  }
  # temp <- stats::rgamma(n=(length(ind)+1), shape=c(m[ind], alpha), scale = 1)
  temp_norm1 <- temp1/sum(temp1)
  temp_norm2 <- temp2/sum(temp2)
  w1[ind1] <- temp_norm1[-length(temp_norm1)]
  R1 <- temp_norm1[length(temp_norm1)]
  w2[ind2] <- temp_norm2[-length(temp_norm2)]
  R2 <- temp_norm2[length(temp_norm2)]
  
  
  #R is the rest, i.e. the weight for potential new clusters
  
  # Sample the latent u
  u1  <- stats::runif(maxCl1)*w1[c1]
  u2  <- stats::runif(maxCl2)*w2[c2]
  u_star1 <- min(u1)
  u_star2 <- min(u2)
  
  # Sample the remaining weights that are needed with stick-breaking
  # i.e. the new clusters
  ind_new <- which(m==0) # potential new clusters
  if(length(ind_new)>0){
    t <- 0 # the number of new non empty clusters
    while(R>u_star && (t<length(ind_new))){
      # sum(w)<1-min(u) <=> R>min(u) car R=1-sum(w)
      t <- t+1
      beta_temp <- stats::rbeta(n=1, shape1=1, shape2=alpha)
      # weight of the new cluster
      w[ind_new[t]] <- R*beta_temp
      R <- R * (1-beta_temp) # remaining weight
    }
    ind_new <- ind_new[1:t]
    
    # Sample the centers and spread of each new cluster from prior
    for (i in 1:t){
      NiW <- rNiW(hyperG0, diagVar)
      U_mu[, ind_new[i]] <- NiW[["mu"]]
      U_Sigma[, , ind_new[i]] <- NiW[["S"]]
    }
  }
  
  fullCl_ind <- which(w != 0)
  
  # likelihood of belonging to each cluster computation
  # sampling clusters
  if(length(fullCl_ind)>1){
    U_mu_full <- sapply(fullCl_ind, function(j) U_mu[, j])
    U_Sigma_list <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
    l <- mmvnpdfC(z, mean=U_mu_full, varcovM=U_Sigma_list, Log = TRUE)
    u_mat <- t(sapply(w[fullCl_ind], function(x){as.numeric(u < x)}))
    prob_mat_log <- log(u_mat) + l
    
    #fast C++ code
    c <- fullCl_ind[sampleClassC(probMat = prob_mat_log, Log = TRUE)]
    #         #slow C++ code
    #         c <- fullCl_ind[sampleClassC_bis(prob_mat)]
    #        #vectorized R code
    #        c <- fullCl_ind[apply(X= prob_mat, MARGIN=2, FUN=function(v){match(1,rmultinom(n=1, size=1, prob=v))})]
    #         #alternative implementation:
    #         prob_colsum <- colSums(prob_mat)
    #         prob_norm <- apply(X=prob_mat, MARGIN=1, FUN=function(r){r/prob_colsum})
    #         c <- fullCl_ind[apply(X=prob_norm, MARGIN=1, FUN=function(r){match(TRUE,stats::runif(1) <cumsum(r))})]
  }else{
    c <- rep(fullCl_ind, maxCl)
  }
  m_new <- numeric(maxCl) # number of observations in each cluster
  m_new[unique(c)] <- table(c)[as.character(unique(c))]
  
  # non vectorized code for cluster allocation:
  #     nb_fullCl <- nb_fullCl + t
  #     l <- numeric(length(fullCl_ind)) # likelihood of belonging to each cluster
  #     m_new <- numeric(maxCl) # number of observations in each cluster
  #     for(i in 1:maxCl){
  #         for (j in fullCl_ind){
  #             l[j] <- mvnpdf(x = matrix(z[,i], ncol= 1, nrow=length(z[,i])) ,
  #                            mean = U_mu[, j],
  #                            varcovM = U_Sigma[, , j])*w[j]
  #         }
  #         c[i] <- rmultinom(n=1, size=1, prob=l)
  #         m_new[c[i]] <- m_new[c[i]] + 1
  #     }
  
  return(list("c"=c, "m"=m_new, "weights"=w,"U_mu"=U_mu,"U_Sigma"=U_Sigma))
}
state1 <- rinit()
state1 <- single_kernel(state1)
state2 <- rinit()

state1$alpha
state2$alpha

##
# nrep <- 10000
# alphas1_ <- rep(0, nrep)
# for (irep in 1:nrep) alphas1_[irep] <- NPflow::sample_alpha(alpha_old=state1$alpha, n=n, K=nbClust1, a=a, b=b)
# alphas2_ <- rep(0, nrep)
# for (irep in 1:nrep) alphas2_[irep] <- NPflow::sample_alpha(alpha_old=state2$alpha, n=n, K=nbClust2, a=a, b=b)
# 
# hist(alphas1_, prob = TRUE, nclass = 60)
# hist(alphas2_, prob = TRUE, nclass = 60, add = TRUE, col = rgb(1,0,0,0.5))
# 
# pairs_ <- matrix(nrow = nrep, ncol = 2)
# for (irep in 1:nrep){
#   pairs_[irep,] <- sample_coupled_alpha(alpha_old1 = state1$alpha, alpha_old2 = state2$alpha, n = n, K1 = nbClust1, K2 = nbClust2, a = a, b = b)
# }
# 
# hist(alphas1_, prob = TRUE, nclass = 60)
# hist(pairs_[,1], prob = TRUE, nclass = 60, add = TRUE, col = rgb(1,0,0,0.5))
# 
# hist(alphas2_, prob = TRUE, nclass = 60)
# hist(pairs_[,2], prob = TRUE, nclass = 60, add = TRUE, col = rgb(1,0,0,0.5))
# 
# mean(apply(pairs_, 1, function(row) row[1]==row[2]))


coupled_kernel <- function(state1, state2){
  nbClust1 <- length(unique(state1$c))
  nbClust2 <- length(unique(state2$c))
  # update alpha
  # alpha1 <- NPflow:::sample_alpha(alpha_old=state1$alpha, n=n, K=nbClust1, a=a, b=b)
  # alpha2 <- NPflow:::sample_alpha(alpha_old=state2$alpha, n=n, K=nbClust2, a=a, b=b)
  alphas <- sample_coupled_alpha(alpha_old1 = state1$alpha, alpha_old2 = state2$alpha, n = n, K1 = nbClust1, K2 = nbClust2, a = a, b = b)
  alpha1 <- alphas[1]
  alpha2 <- alphas[2]
  # slice sampling update (c)
  slice1 <- NPflow:::sliceSampler_N(c=state1$c, m=state1$m, alpha=state1$alpha, z=z, hyperG0=hyperG0,
                                   U_mu=state1$U_mu, U_Sigma=state1$U_Sigma, diagVar=diagVar)
  slice2 <- NPflow:::sliceSampler_N(c=state2$c, m=state2$m, alpha=state2$alpha, z=z, hyperG0=hyperG0,
                                   U_mu=state2$U_mu, U_Sigma=state2$U_Sigma, diagVar=diagVar)
  m <- slice[["m"]]
  c <- slice[["c"]]
  # weights_list[[i]] <- slice[["weights"]]
  U_mu<-slice[["U_mu"]]
  U_Sigma<-slice[["U_Sigma"]]
  # Update cluster locations
  fullCl <- which(m!=0)
  for(j in fullCl){
    obs_j <- which(c==j)
    #cat("cluster ", j, ":\n")
    if (use_variance_hyperprior){
      U_SS[[j]] <- NPflow:::update_SS(z=z[, obs_j, drop=FALSE], S=hyperG0, hyperprior = list("Sigma"=U_Sigma[,,j]))
    } else{
      U_SS[[j]] <- NPflow:::update_SS(z=z[, obs_j, drop=FALSE], S=hyperG0)
    }
    NiW <- NPflow:::rNiW(U_SS[[j]], diagVar=diagVar)
    U_mu[, j] <- NiW[["mu"]]
    U_SS[[j]][["mu"]] <- NiW[["mu"]]
    
    U_Sigma[, , j] <- NiW[["S"]]
    U_SS[[j]][["S"]] <- NiW[["S"]]
    
    # U_SS[[j]][["weight"]] <- weights_list[[i]][j]
    #cat("sampled S =", NiW[["S"]], "\n\n\n")
  }
  return(list(alpha = alpha, c = c, U_mu = U_mu, U_Sigma = U_Sigma, U_SS = U_SS))
}


