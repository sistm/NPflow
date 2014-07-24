#'Convergence diagnostic plots
#'
#'
#'@export

plot_ConvDPM <- function(MCMCsample, from=1, to=length(MCMCsample$logposterior_list), shift=0){
    
    
    from_lab <- from + shift
    to_lab <- to + shift
    
    par("mfrow"=c(3,2))
    
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 1))[from:to],
         x=c(from_lab:to_lab), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior data")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 2))[from:to], 
         x=c(from_lab:to_lab), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior partition")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 3))[from:to], 
         x=c(from_lab:to_lab), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior NiW prior")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 4))[from:to], 
         x=c(from_lab:to_lab), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior alpha prior")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, 'sum'))[from:to], 
         x=c(from_lab:to_lab), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior All")
}