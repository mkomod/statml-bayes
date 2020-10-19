library(Rcpp)
library(mvtnorm)

Rcpp::sourceCpp("../rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

P <- 1
theta.p <- c(1)
NS <- c(1e2, 1e3, 1e4, 2.5e4)

res <- data.frame(N=numeric(), par=numeric(), val=numeric(), mse=numeric())

for (i in seq_along(NS)) {
    N <- NS[i]
    theta.init <- c(2)
    X <- matrix(rlaplace(N * P), ncol=P)
    Y <- matrix(rmvnorm(N, mean=rep(0, P)), ncol=P)
    opt <- optim(theta.init, fn=J, X=X, Y=Y, method="CG",
		 control=list(maxit=1e3))
    res <- rbind(res, 
		 list(N=N, par=opt$par, val=opt$val, 
		      mse=mean(theta.p - opt$par)^2))
}

plot(log10(res$N), log10(res$mse))

