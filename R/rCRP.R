#'Generating cluster data from the Chinese Restaurant Process
#'
#'@param n
#'
#'@param alpha
#'
#'@param hyperG0
#'
#'@export rCRP
#'
#'@examples
#'f <- function(x){sum(x/(x+0:12751))-39}
#'alpha_otpim <- uniroot(f, c(0.0000001,100))
#'alpha_otpim$root
#'set.seed(50)
#'GvHD_sims <- rCRP(n=12751, alpha=alpha_otpim$root, hyperG0)
#'s <- summary(factor(GvHD_sims$cluster))
#'s_ord <- s[order(s, decreasing = TRUE)]
#'names(s_ord)<-1:length(s_ord)
#'barplot(s_ord,
#'xlab="Cluster", ylab="Size", ylim=c(0,5000),
#'main="DP(alpha=4.9): CRP simulations\n(n=12,751)")
#'
#'plot(y=as.numeric(names(summary(factor(s_ord)))), 
#'x=summary(factor(s_ord)), 
#'log="xy", col="blue", pch=20, ylim=c(1,5000),
#'main="DP(alpha=4.9): CRP simulations\n(n=12,751)", 
#'xlab="Number of clusters", ylab="Cluster size")
#'
#'s_hist <- hist(log(s), col="grey")
#'plot(x=exp(s_hist$mids), y=s_hist$counts, pch=20,log="xy", 
#'     xlab="Cluster size", ylab="Number of clusters",col="blue", type="p", cex=1.5)
#'
#'barplot(height=s_hist$count,log="y", 
#'        xlab="Cluster size", ylab="Number of clusters", space=0,
#'        main="DP(alpha=4.9): CRP simulations\n(n=12,751)")
#'axis(1, at= log(c(2,5,10,20,50,100,200,500,1000,2000,5000)), 
#'     labels=c(2,5,10,20,50,100,200,500,1000,2000,5000))
#'
rCRP <- function(n=1000, alpha=10, hyperG0){
    
    theta <- list()
    cluster <- numeric(n)
        
    for (c in 1:n){
        p0 <- alpha/(c-1+alpha)
        u <- runif(n=1, min = 0, max = 1)
        
        if (u<p0){
            # Accept: sample new value
            # cat("acceptation:", u, "<", p0, "\n")
            cluster[c] <- max(cluster)+1
            theta[[cluster[c]]] <- nniw_rnd(hyperG0, diagVar=FALSE)
            
        }else{
            # Reject: sample old value
            # cat("rejection:", u, ">=", p0, "\n")
            u1  <-  u - p0
            weights <- summary(factor(cluster))[-1]
            cluster[c] <- which(rmultinom(n=1, size=1, prob=weights)==1)
        }
        cat(c,"/", n," sim\n", sep="")
    }
    
    return(list("theta"=theta, "cluster"=cluster))
}