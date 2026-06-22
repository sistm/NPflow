# Generating cluster data from the Chinese Restaurant Process

Generating cluster data from the Chinese Restaurant Process

## Usage

``` r
rCRP(n = 1000, alpha = 2, hyperG0, verbose = TRUE)
```

## Arguments

- n:

  number of observations

- alpha:

  concentration parameter

- hyperG0:

  base distribution hyperparameter

- verbose:

  logical flag indicating whether info is written in the console.

## Examples

``` r

rm(list=ls())

d=2
hyperG0 <- list()
hyperG0[["NNiW"]] <- list()
hyperG0[["NNiW"]][["b_xi"]] <- rep(0,d)
hyperG0[["NNiW"]][["b_psi"]] <- rep(0,d)
hyperG0[["NNiW"]][["D_xi"]] <- 100
hyperG0[["NNiW"]][["D_psi"]] <- 8
hyperG0[["NNiW"]][["nu"]] <- d+1
hyperG0[["NNiW"]][["lambda"]] <- diag(c(1,1))

hyperG0[["scale"]] <- list()

set.seed(4321)
N <- 200
alph <- runif(n=1,0.2,2)
GvHD_sims <- rCRP(n=2*N, alpha=alph, hyperG0=hyperG0)
#> 1/400 sim
#> 2/400 sim
#> 3/400 sim
#> 4/400 sim
#> 5/400 sim
#> 6/400 sim
#> 7/400 sim
#> 8/400 sim
#> 9/400 sim
#> 10/400 sim
#> 11/400 sim
#> 12/400 sim
#> 13/400 sim
#> 14/400 sim
#> 15/400 sim
#> 16/400 sim
#> 17/400 sim
#> 18/400 sim
#> 19/400 sim
#> 20/400 sim
#> 21/400 sim
#> 22/400 sim
#> 23/400 sim
#> 24/400 sim
#> 25/400 sim
#> 26/400 sim
#> 27/400 sim
#> 28/400 sim
#> 29/400 sim
#> 30/400 sim
#> 31/400 sim
#> 32/400 sim
#> 33/400 sim
#> 34/400 sim
#> 35/400 sim
#> 36/400 sim
#> 37/400 sim
#> 38/400 sim
#> 39/400 sim
#> 40/400 sim
#> 41/400 sim
#> 42/400 sim
#> 43/400 sim
#> 44/400 sim
#> 45/400 sim
#> 46/400 sim
#> 47/400 sim
#> 48/400 sim
#> 49/400 sim
#> 50/400 sim
#> 51/400 sim
#> 52/400 sim
#> 53/400 sim
#> 54/400 sim
#> 55/400 sim
#> 56/400 sim
#> 57/400 sim
#> 58/400 sim
#> 59/400 sim
#> 60/400 sim
#> 61/400 sim
#> 62/400 sim
#> 63/400 sim
#> 64/400 sim
#> 65/400 sim
#> 66/400 sim
#> 67/400 sim
#> 68/400 sim
#> 69/400 sim
#> 70/400 sim
#> 71/400 sim
#> 72/400 sim
#> 73/400 sim
#> 74/400 sim
#> 75/400 sim
#> 76/400 sim
#> 77/400 sim
#> 78/400 sim
#> 79/400 sim
#> 80/400 sim
#> 81/400 sim
#> 82/400 sim
#> 83/400 sim
#> 84/400 sim
#> 85/400 sim
#> 86/400 sim
#> 87/400 sim
#> 88/400 sim
#> 89/400 sim
#> 90/400 sim
#> 91/400 sim
#> 92/400 sim
#> 93/400 sim
#> 94/400 sim
#> 95/400 sim
#> 96/400 sim
#> 97/400 sim
#> 98/400 sim
#> 99/400 sim
#> 100/400 sim
#> 101/400 sim
#> 102/400 sim
#> 103/400 sim
#> 104/400 sim
#> 105/400 sim
#> 106/400 sim
#> 107/400 sim
#> 108/400 sim
#> 109/400 sim
#> 110/400 sim
#> 111/400 sim
#> 112/400 sim
#> 113/400 sim
#> 114/400 sim
#> 115/400 sim
#> 116/400 sim
#> 117/400 sim
#> 118/400 sim
#> 119/400 sim
#> 120/400 sim
#> 121/400 sim
#> 122/400 sim
#> 123/400 sim
#> 124/400 sim
#> 125/400 sim
#> 126/400 sim
#> 127/400 sim
#> 128/400 sim
#> 129/400 sim
#> 130/400 sim
#> 131/400 sim
#> 132/400 sim
#> 133/400 sim
#> 134/400 sim
#> 135/400 sim
#> 136/400 sim
#> 137/400 sim
#> 138/400 sim
#> 139/400 sim
#> 140/400 sim
#> 141/400 sim
#> 142/400 sim
#> 143/400 sim
#> 144/400 sim
#> 145/400 sim
#> 146/400 sim
#> 147/400 sim
#> 148/400 sim
#> 149/400 sim
#> 150/400 sim
#> 151/400 sim
#> 152/400 sim
#> 153/400 sim
#> 154/400 sim
#> 155/400 sim
#> 156/400 sim
#> 157/400 sim
#> 158/400 sim
#> 159/400 sim
#> 160/400 sim
#> 161/400 sim
#> 162/400 sim
#> 163/400 sim
#> 164/400 sim
#> 165/400 sim
#> 166/400 sim
#> 167/400 sim
#> 168/400 sim
#> 169/400 sim
#> 170/400 sim
#> 171/400 sim
#> 172/400 sim
#> 173/400 sim
#> 174/400 sim
#> 175/400 sim
#> 176/400 sim
#> 177/400 sim
#> 178/400 sim
#> 179/400 sim
#> 180/400 sim
#> 181/400 sim
#> 182/400 sim
#> 183/400 sim
#> 184/400 sim
#> 185/400 sim
#> 186/400 sim
#> 187/400 sim
#> 188/400 sim
#> 189/400 sim
#> 190/400 sim
#> 191/400 sim
#> 192/400 sim
#> 193/400 sim
#> 194/400 sim
#> 195/400 sim
#> 196/400 sim
#> 197/400 sim
#> 198/400 sim
#> 199/400 sim
#> 200/400 sim
#> 201/400 sim
#> 202/400 sim
#> 203/400 sim
#> 204/400 sim
#> 205/400 sim
#> 206/400 sim
#> 207/400 sim
#> 208/400 sim
#> 209/400 sim
#> 210/400 sim
#> 211/400 sim
#> 212/400 sim
#> 213/400 sim
#> 214/400 sim
#> 215/400 sim
#> 216/400 sim
#> 217/400 sim
#> 218/400 sim
#> 219/400 sim
#> 220/400 sim
#> 221/400 sim
#> 222/400 sim
#> 223/400 sim
#> 224/400 sim
#> 225/400 sim
#> 226/400 sim
#> 227/400 sim
#> 228/400 sim
#> 229/400 sim
#> 230/400 sim
#> 231/400 sim
#> 232/400 sim
#> 233/400 sim
#> 234/400 sim
#> 235/400 sim
#> 236/400 sim
#> 237/400 sim
#> 238/400 sim
#> 239/400 sim
#> 240/400 sim
#> 241/400 sim
#> 242/400 sim
#> 243/400 sim
#> 244/400 sim
#> 245/400 sim
#> 246/400 sim
#> 247/400 sim
#> 248/400 sim
#> 249/400 sim
#> 250/400 sim
#> 251/400 sim
#> 252/400 sim
#> 253/400 sim
#> 254/400 sim
#> 255/400 sim
#> 256/400 sim
#> 257/400 sim
#> 258/400 sim
#> 259/400 sim
#> 260/400 sim
#> 261/400 sim
#> 262/400 sim
#> 263/400 sim
#> 264/400 sim
#> 265/400 sim
#> 266/400 sim
#> 267/400 sim
#> 268/400 sim
#> 269/400 sim
#> 270/400 sim
#> 271/400 sim
#> 272/400 sim
#> 273/400 sim
#> 274/400 sim
#> 275/400 sim
#> 276/400 sim
#> 277/400 sim
#> 278/400 sim
#> 279/400 sim
#> 280/400 sim
#> 281/400 sim
#> 282/400 sim
#> 283/400 sim
#> 284/400 sim
#> 285/400 sim
#> 286/400 sim
#> 287/400 sim
#> 288/400 sim
#> 289/400 sim
#> 290/400 sim
#> 291/400 sim
#> 292/400 sim
#> 293/400 sim
#> 294/400 sim
#> 295/400 sim
#> 296/400 sim
#> 297/400 sim
#> 298/400 sim
#> 299/400 sim
#> 300/400 sim
#> 301/400 sim
#> 302/400 sim
#> 303/400 sim
#> 304/400 sim
#> 305/400 sim
#> 306/400 sim
#> 307/400 sim
#> 308/400 sim
#> 309/400 sim
#> 310/400 sim
#> 311/400 sim
#> 312/400 sim
#> 313/400 sim
#> 314/400 sim
#> 315/400 sim
#> 316/400 sim
#> 317/400 sim
#> 318/400 sim
#> 319/400 sim
#> 320/400 sim
#> 321/400 sim
#> 322/400 sim
#> 323/400 sim
#> 324/400 sim
#> 325/400 sim
#> 326/400 sim
#> 327/400 sim
#> 328/400 sim
#> 329/400 sim
#> 330/400 sim
#> 331/400 sim
#> 332/400 sim
#> 333/400 sim
#> 334/400 sim
#> 335/400 sim
#> 336/400 sim
#> 337/400 sim
#> 338/400 sim
#> 339/400 sim
#> 340/400 sim
#> 341/400 sim
#> 342/400 sim
#> 343/400 sim
#> 344/400 sim
#> 345/400 sim
#> 346/400 sim
#> 347/400 sim
#> 348/400 sim
#> 349/400 sim
#> 350/400 sim
#> 351/400 sim
#> 352/400 sim
#> 353/400 sim
#> 354/400 sim
#> 355/400 sim
#> 356/400 sim
#> 357/400 sim
#> 358/400 sim
#> 359/400 sim
#> 360/400 sim
#> 361/400 sim
#> 362/400 sim
#> 363/400 sim
#> 364/400 sim
#> 365/400 sim
#> 366/400 sim
#> 367/400 sim
#> 368/400 sim
#> 369/400 sim
#> 370/400 sim
#> 371/400 sim
#> 372/400 sim
#> 373/400 sim
#> 374/400 sim
#> 375/400 sim
#> 376/400 sim
#> 377/400 sim
#> 378/400 sim
#> 379/400 sim
#> 380/400 sim
#> 381/400 sim
#> 382/400 sim
#> 383/400 sim
#> 384/400 sim
#> 385/400 sim
#> 386/400 sim
#> 387/400 sim
#> 388/400 sim
#> 389/400 sim
#> 390/400 sim
#> 391/400 sim
#> 392/400 sim
#> 393/400 sim
#> 394/400 sim
#> 395/400 sim
#> 396/400 sim
#> 397/400 sim
#> 398/400 sim
#> 399/400 sim
#> 400/400 sim
library(ggplot2)
q <- (ggplot(data=cbind.data.frame("D1"=GvHD_sims$data[1,],
                                  "D2"=GvHD_sims$data[2,],
                                  "Cluster"=GvHD_sims$cluster),
             aes(x=D1, y=D2))
      + geom_point(aes(colour=Cluster), alpha=0.6)
      + theme_bw()
      )
q

#q + stat_density2d(alpha=0.15, geom="polygon")

if(interactive()){
MCMCy1 <- DPMGibbsSkewT(z=GvHD_sims$data[,1:N],
                        hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
                        doPlot=TRUE, nbclust_init=64, plotevery=500,
                        gg.add=list(theme_bw()), diagVar=FALSE)
 s1 <- summary(MCMCy1, burnin=4000, thin=5,
               posterior_approx=TRUE)
 F1 <- FmeasureC(ref=GvHD_sims$cluster[1:N], pred=s1$point_estim$c_est)

 # s <- summary(MCMCy1, burnin=4000, thin=5,
 #               posterior_approx=TRUE, K=1)
 # s2 <- summary(MCMCy1, burnin=4000, thin=5,
 #               posterior_approx=TRUE, K=2)
 # MCMCy2_seqPost<- DPMGibbsSkewT(z=GvHD_sims$data[,(N+1):(2*N)],
 #                                  hyperG0=s1$param_post$parameters,
 #                                  a=s1$param_post$alpha_param$shape,
 #                                  b=s1$param_post$alpha_param$rate,
 #                                  N=5000, doPlot=TRUE, nbclust_init=64, plotevery=500,
 #                                  gg.add=list(theme_bw()), diagVar=FALSE)

 MCMCy2_seqPost <- DPMGibbsSkewT_SeqPrior(z=GvHD_sims$data[,(N+1):(2*N)],
                                           prior=s1$param_post, hyperG0=hyperG0$NNiW, , N=1000,
                                           doPlot=TRUE, nbclust_init=10, plotevery=100,
                                           gg.add=list(theme_bw()), diagVar=FALSE)
 s2_seqPost <- summary(MCMCy2_seqPost, burnin=600, thin=2)
 F2_seqPost <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=s2_seqPost$point_estim$c_est)

 MCMCy2 <- DPMGibbsSkewT(z=GvHD_sims$data[,(N+1):(2*N)],
                         hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
                         doPlot=TRUE, nbclust_init=64, plotevery=500,
                         gg.add=list(theme_bw()), diagVar=FALSE)
 s2 <- summary(MCMCy2, burnin=4000, thin=5)
 F2 <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=s2$point_estim$c_est)

 MCMCtot <- DPMGibbsSkewT(z=GvHD_sims$data,
                          hyperG0$NNiW, a=0.0001, b=0.0001, N=5000,
                          doPlot=TRUE, nbclust_init=10, plotevery=500,
                          gg.add=list(theme_bw()), diagVar=FALSE)
 stot <- summary(MCMCtot, burnin=4000, thin=5)
 F2tot <- FmeasureC(ref=GvHD_sims$cluster[(N+1):(2*N)], pred=stot$point_estim$c_est[(N+1):(2*N)])

 c(F1, F2, F2_seqPost, F2tot)
}
```
