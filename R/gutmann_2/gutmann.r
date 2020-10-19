library(Rcpp)
library(mvtnorm)

Rcpp::sourceCpp("../rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

NS <- round(10^seq(2, 4.5, length.out=15))
res <- data.frame(N=numeric(), par=numeric(), val=numeric(), mse=numeric())

P <- 1
theta.p <- c(1/sqrt(2))
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

plot(log10(res$N), log10(res$mse), log='x')

P <- 2
theta.p <- c(0.5)
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


P <- 3
theta.p <- c(1/(sqrt(2)^P))
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

