#'Convergence diagnostic plots
#'
#'
#'@param thin integer giving the spacing at which MCMC iterations are kept.
#'Default is \code{1}, i.e. no thining.
#'
#'@param ... further arguments passed to or from other methods
#'
#'@export

plot_ConvDPM <- function(MCMCsample, from=1, to=length(MCMCsample$logposterior_list),
                         shift=0, thin=1, ...){

    from_lab <- from + shift
    to_lab <- to*thin + shift

    par("mfrow"=c(3,2))

    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 1))[from:to],
         x=seq(from=from_lab, to=to_lab, by=thin), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior data", ...)
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 2))[from:to],
         x=seq(from=from_lab, to=to_lab, by=thin), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior partition",...)
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 3))[from:to],
         x=seq(from=from_lab, to=to_lab, by=thin), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior NiW prior",...)
    plot(y=unlist(lapply(MCMCsample$logposterior_list, '[', 4))[from:to],
         x=seq(from=from_lab, to=to_lab, by=thin), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior alpha prior",...)
    plot(y=unlist(lapply(MCMCsample$logposterior_list, 'sum'))[from:to],
         x=seq(from=from_lab, to=to_lab, by=thin), type="b", col="blue", pch=8, cex=0.7,
         xlab="MCMC iteration", ylab="log posterior All",...)
}