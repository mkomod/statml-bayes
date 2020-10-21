library(Rcpp)
library(mvtnorm)
library(parallel)

CORES <- 4

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

P <- 4
NS <- round(10^seq(2, 2.5, length.out=CORES))
res <- mclapply(NS, function(N) {
    r <- data.frame(N=numeric(), run=numeric(), s=numeric(), val=numeric(), 
		    mse=numeric(), t=numeric())
    set.seed(123)
    for (run in 1:200) {
	A <- matrix(runif(16, 0, 10), nrow=4)
	S <- A %*% t(A)
	X <- matrix(rlaplace(N * P), ncol=P) %*% A
	Y <- matrix(rmvnorm(N, mean=rep(0, P), S), ncol=P)
	thetaP <- c(as.vector(matrix(solve(A),ncol=1)), 
		    log(abs(det(solve(A)))) - log(4))

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
    cat(N, "done\n")
    # assign(paste0("r_", N), r)
    # save(list=paste0("r_", N), file=paste0("r_", N, ".RData"))
    return(r)
}, mc.cores=length(NS))

# MSE
m <- sapply(res, function(r) r$mse)
mm <- apply(m, 2, mean)
plot(log10(mm))
polygon(c(NS, rev(NS)), log10(c(res$upr, rev(res$lwr))), 
	border=NA, col=rgb(0, 1, 0, 0.3))

# Comp time
tm <- sapply(res, function(r) r$t)
tmm <- apply(tm, 2, mean)
plot(log10(mm), log10(tmm))



