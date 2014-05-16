#'
#'plot_ConvDPM
#'
#'@export plot_ConvDPM

plot_ConvDPM <- function(MCMCsample, from=1, to=length(MCMCsample$logposterior_list)){
    
    par("mfrow"=c(3,2))
    
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 1))[from:to],
         x=c(from:to), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior data")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 2))[from:to], 
         x=c(from:to), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior partition")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 3))[from:to], 
         x=c(from:to), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior NiW prior")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 4))[from:to], 
         x=c(from:to), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior alpha prior")
    plot(y=unlist(lapply(MCMCsample$logposterior_list, 'sum'))[from:to], 
         x=c(from:to), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior All")
}