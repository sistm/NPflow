# # Stick-breaking realizations
# strickBreak_p = function(v) {
#     cumv = cumprod(1 - v)
#     
#     cumv_shift <- c(1, cumv[-length(v)])
#     p <- v*cumv_shift
#     
# #     p = numeric(length(v))
# #     p[1] = v[1]
# #     for (i in 2:length(v)){
# #         p[i] = v[i] * cumv[i - 1]
# #     }
#     return(p)
# }
# 
# rDP = function(alpha=10, G0=rnorm, N = 1000) {
#     theta = G0(N)
#     v = rbetInt(N, 1, alpha)
#     p = strickBreak_p(v)
#     return(list("G0"=G0, "theta" = theta, "v" = v, "p" = p))
# }
# 
# plot_rDP = function(rdp, ...) {
#     plot(rdp$theta, rdp$p, type = "n", ylim=c(0, 0.5),
#          xlab="x", ylab="Density", ...)
#     curve(dnorm, col="blue", lwd=3, add=T)
#     lines(rdp$theta, rdp$p, type = "h", ylim=c(0, 0.5), lwd=3)
#     lines(rdp$theta, rdp$p, type = "p", pch=16)
# }
# 
# plot_rDP(rDP(alpha=3))
# 
# 
# 
# 
# # Examples plots ----
# set.seed(15)
# rdp <- rDP(alpha=2)
# pdf(width=10, height=8, file="ExDPM_alpha2.pdf")
#     plot_rDP(rdp, main="alpha=2", cex.main=1.5)
#     e<- NULL
#     for(i in 1:length(rdp$p)){
#         nb <- floor(500*rdp$p[i])
#         if(nb>0){
#             e <- c(e, rnorm(n=nb, m=rdp$theta[i]), sd=0.01)
#         }
#     }
#     lines(density(e, bw=0.2), lty=2, lwd=3)
#     legend("topleft", 
#            c("Sample G~DP(alpha=2, G0)",
#              "Base distribution G0=N(0,1)", 
#              "Dirichlet Process Mixture of Gaussians"),
#            lwd=2, pch=16, lty=c(0,1,2), col=c("black", "blue", "black"),
#            cex=1.1, pt.cex=c(1,0,0))
# dev.off()
# 
# set.seed(15)
# rdp <- rDP(alpha=20)
# pdf(width=10, height=8, file="ExDPM_alpha20.pdf")
#     plot_rDP(rdp, main="alpha=20", cex.main=1.5)
#     e<- NULL
#     for(i in 1:length(rdp$p)){
#         nb <- floor(500*rdp$p[i])
#         if(nb>0){
#             e <- c(e, rnorm(n=nb, m=rdp$theta[i]), sd=0.01)
#         }
#     }
#     lines(density(e, bw=0.2), lty=2, lwd=3)
#     legend("topleft", 
#            c("Sample G~DP(alpha=2, G0)",
#              "Base distribution G0=N(0,1)", 
#              "Dirichlet Process Mixture of Gaussians"),
#            lwd=2, pch=16, lty=c(0,1,2), col=c("black", "blue", "black"),
#            cex=1.1, pt.cex=c(1,0,0))
# dev.off()
# 
# 
# 
# 
# # Slice sampler ----
# set.seed(15)
# rdp <- rDP(alpha=5)
# cut1 <- floor(100*rdp$p[order(rdp$p, decreasing=T)][5])/100
# cut2 <- floor(100*rdp$p[order(rdp$p, decreasing=T)][9])/100
# pdf(width=8, height=6.5, file="SliceDP_1.pdf")
#     plot(rdp$theta, rdp$p, type = "n", xlab="theta", ylab="DPM density")
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           type = "h", lwd=2)
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           pch=16, type = "p", cex=0.9)
# dev.off()
# 
# pdf(width=8, height=6.5, file="SliceDP_2.pdf")
#     plot(rdp$theta, rdp$p, type = "n", xlab="theta", ylab="DPM density")
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           type = "h", lwd=2, col="blue")
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           pch=16, type = "p", cex=0.9, col="blue")
#     lines(rdp$theta[which(rdp$p<cut1)], rdp$p[which(rdp$p<cut1)], 
#           type = "h", lwd=1, col="black")
#     lines(rdp$theta[which(rdp$p<cut1)], rdp$p[which(rdp$p<cut1)], 
#           pch=16, type = "p", cex=0.6)
# dev.off()
# 
# pdf(width=8, height=6.5, file="SliceDP_3.pdf")
#     plot(rdp$theta, rdp$p, type = "n", xlab="theta", ylab="DPM density")
#     abline(col="grey70", h=cut2, lwd=3)
#     text(paste("u =", cut2), x=-4.3, y=0.05, col="grey70", cex=1.4, adj=0)
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           type = "h", lwd=2, col="blue")
#     lines(rdp$theta[which(rdp$p>cut1)], rdp$p[which(rdp$p>cut1)], 
#           pch=16, type = "p", cex=0.9, col="blue")
# lines(rdp$theta[which(rdp$p>cut2 & rdp$p<cut1)], rdp$p[which(rdp$p>cut2 & rdp$p<cut1)], 
#       type = "h", lwd=2, col="black")
# lines(rdp$theta[which(rdp$p>cut2 & rdp$p<cut1)], rdp$p[which(rdp$p>cut2 & rdp$p<cut1)], 
#       pch=16, type = "p", cex=0.9, col="black")
#     lines(rdp$theta[which(rdp$p<cut2)], rdp$p[which(rdp$p<cut2)], 
#           type = "h", lwd=1, col="grey")
#     lines(rdp$theta[which(rdp$p<cut2 )], rdp$p[which(rdp$p<cut2)], 
#           pch=16, type = "p", cex=0.6, col="grey")
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# # Slice sampler ----
# slice = function(n, init_x, target, Int) {
#     u = x = rep(NA, n)
#     x[1] = init_x
#     u[1] = runif(1, 0, target(x[1]))  # This never actually gets used
#     
#     for (i in 2:n) {
#         u[i] = runif(1, 0, target(x[i - 1]))
#         endpoints = Int(u[i], x[i - 1])  # The second argument is used in the second example
#         x[i] = runif(1, endpoints[1], endpoints[2])
#     }
#     return(list(x = x, u = u))
# }
# 
# # Normal slice sampling
# Int = function(u, xx) {
#     c(uniroot(function(x) dnorm(x) - u, c(-10^10, xx))$root, uniroot(function(x) dnorm(x) - 
#                                                                          u, c(xx, 10^10))$root)
# }
# 
# set.seed(6)
# res = slice(10, 0.1, dnorm, A)
# x = res$x
# u = res$u
# 
# i = 1
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=0", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="blue")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=1", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="blue")
# segments(x[i], 0, x[i], dnorm(x[i]), col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=1", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="blue")
# segments(x[i], 0, x[i], dnorm(x[i]), col = "gray")
# arrows(x[i], u[i], x[i], u[i + 1], length = 0.1)
# points(x[i], u[i + 1])
# segments(Int(u[i + 1], x[i])[1], u[i + 1], Int(u[i + 1], x[i])[2], u[i + 1], col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=1", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# segments(Int(u[i + 1], x[i])[1], u[i + 1], Int(u[i + 1], x[i])[2], u[i + 1], col = "gray")
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="blue")
# 
# j = 2
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=2", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="blue")
# segments(x[j], 0, x[j], dnorm(x[j]), col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=2", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="blue")
# segments(x[j], 0, x[j], dnorm(x[j]), col = "gray")
# arrows(x[j], u[j], x[j], u[j + 1], length = 0.1)
# points(x[j], u[j + 1])
# segments(Int(u[j + 1], x[j])[1], u[j + 1], Int(u[j + 1], x[j])[2], u[j + 1], col = "gray")
# 
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=2", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# segments(Int(u[j + 1], x[j])[1], u[j + 1], Int(u[j + 1], x[j])[2], u[j + 1], col = "gray")
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="blue")
# 
# 
# k=3
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=3", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="blue")
# segments(x[k], 0, x[k], dnorm(x[k]), col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=3", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="blue")
# segments(x[k], 0, x[k], dnorm(x[k]), col = "gray")
# arrows(x[k], u[k], x[k], u[k + 1], length = 0.1)
# points(x[k], u[k + 1])
# segments(Int(u[k + 1], x[k])[1], u[k + 1], Int(u[k + 1], x[k])[2], u[k + 1], col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=3", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="black")
# segments(x[k], u[k], x[k], u[k + 1], length = 0.1)
# segments(Int(u[k + 1], x[k])[1], u[k + 1], Int(u[k + 1], x[k])[2], u[k + 1], col = "gray")
# arrows(x[k], u[k + 1], x[k + 1], u[k + 1], length = 0.1)
# points(x[k + 1], u[k + 1], pch = 19, col="blue")
# 
# l=4
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=4", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="black")
# segments(x[k], u[k], x[k], u[k + 1], length = 0.1)
# arrows(x[k], u[k + 1], x[k + 1], u[k + 1], length = 0.1)
# points(x[k + 1], u[k + 1], pch = 19, col="blue")
# segments(x[l], 0, x[l], dnorm(x[l]), col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=4", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="black")
# segments(x[k], u[k], x[k], u[k + 1], length = 0.1)
# arrows(x[k], u[k + 1], x[k + 1], u[k + 1], length = 0.1)
# points(x[k + 1], u[k + 1], pch = 19, col="blue")
# segments(x[l], 0, x[l], dnorm(x[l]), col = "gray")
# arrows(x[l], u[l], x[l], u[l + 1], length = 0.1)
# points(x[l], u[l + 1])
# segments(Int(u[l + 1], x[l])[1], u[l + 1], Int(u[l + 1], x[l])[2], u[l + 1], col = "gray")
# 
# curve(dnorm, -3, 3, ylim = c(0, 0.5), main="i=4", main.cex=1.5)
# points(x[i], u[i], pch = 19, col="black")
# segments(x[i], u[i], x[i], u[i + 1], length = 0.1)
# arrows(x[i], u[i + 1], x[i + 1], u[i + 1], length = 0.1)
# points(x[i + 1], u[i + 1], pch = 19, col="black")
# segments(x[j], u[j], x[j], u[j + 1], length = 0.1)
# arrows(x[j], u[j + 1], x[j + 1], u[j + 1], length = 0.1)
# points(x[j + 1], u[j + 1], pch = 19, col="black")
# segments(x[k], u[k], x[k], u[k + 1], length = 0.1)
# arrows(x[k], u[k + 1], x[k + 1], u[k + 1], length = 0.1)
# points(x[k + 1], u[k + 1], pch = 19, col="black")
# segments(x[l], u[l], x[l], u[l + 1], length = 0.1)
# segments(Int(u[l + 1], x[l])[1], u[l + 1], Int(u[l + 1], x[l])[2], u[l + 1], col = "gray")
# arrows(x[l], u[l + 1], x[l + 1], u[l + 1], length = 0.1)
# points(x[l + 1], u[l + 1], pch = 19, col="blue")
# 
# 
# 
