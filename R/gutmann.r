library(Rcpp)
library(mvtnorm)
library(RColorBrewer)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

P <- 4
NS <- round(10^seq(2, 4.5, length.out=11))
res <- mclapply(NS, function(N) {
    r <- data.frame(N=numeric(), run=numeric(), s=numeric(), val=numeric(), 
		    mse=numeric(), t=numeric())
    for (run in 1:500) {
	A <- matrix(runif(16, 0, 10), nrow=4)
	S <- A %*% t(A)
	X <- matrix(rlaplace(N * P), ncol=P) %*% A
	Y <- matrix(rmvnorm(N, mean=rep(0, P), S), ncol=P)
	thetaP <- c(t(solve(A)), log(abs(det(solve(A)))) - log(4))

	for (i in 1:5) {
	    A_init = matrix(runif(16, 0,10), nrow = 4)
	    B_init = solve(A)
	    c_init = log(abs(det(B_init))) - log(4)
	    theta_init = c(as.vector(matrix(B_init,ncol=1)),c_init)
	    s <- system.time({ 
		opt <- optim(theta_init, J, method="CG", X=X, Y=Y, S=S)
	    })
	    r <- rbind(r, list(N=N, run=run, s=i, val=opt$val, 
		mse=mean((thetaP - opt$par)^2), t=as.numeric(s)[2]))

	}
    }
    return(r)
}, mc.cores=length(NS))


save(res, file="res.RData")

