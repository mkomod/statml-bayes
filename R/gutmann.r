library(Rcpp)
library(mvtnorm)
library(RColorBrewer)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

NS <- round(10^seq(2, 4.5, length.out=11))
PS <- 1:8

res <- list()
for (P in PS) {
    r <- data.frame(N=numeric(), par=numeric(), 
		      val=numeric(), mse=numeric(), t=numeric())
    theta.p <- c(1/sqrt(2)^P)
    for (N in NS) {
	X <- matrix(rlaplace(N * P), ncol=P)
	Y <- matrix(rmvnorm(N, mean=rep(0, P)), ncol=P)
	for (i in 1:5) {
	    theta.init <- runif(1, min=0, max=5)
	    s <- system.time({ 
	    opt <- optim(theta.init, fn=J, X=X, Y=Y, method="CG",
		control=list(maxit=1e3))
	    })
	    r <- rbind(r, list(N=N, par=opt$par, val=opt$val, 
		mse=mean(theta.p - opt$par)^2, t=as.numeric(s)[2]))
	}
    }
    res[[P]] = r
}

save(res, file="res.RData")


# Plots for dists -------------------------------------------------------------

x <- y <- seq(-3, 3, length.out=600)
nd <- length(x)

d.norm <- matrix(apply(expand.grid(x, y), 1, 
		function(x) pn(c(x[1], x[2]))), nd)
d.lap <- matrix(apply(expand.grid(x, y), 1, 
		function(x) pm(c(x[1], x[2]), 1/sqrt(2))), nd)

col.pal <- colorRampPalette(rev(brewer.pal(9, "Blues")))
cols <- col.pal(35)

filled.contour(x, y, d.norm, xlab=expression(X[1]), 
	       ylab=expression(X[2]), col=cols)
filled.contour(x, y, d.lap, xlab=expression(X[1]), 
	       ylab=expression(X[2]), col=cols)

