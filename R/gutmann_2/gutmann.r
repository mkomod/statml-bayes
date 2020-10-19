library(Rcpp)
library(mvtnorm)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("J.cpp")

P <- 1
theta.p <- c(1)
for (i in c(1e2, 1e3, 1e4, 2.5e4, 5e4)) {
    N <- i
    theta.init <- c(2)
    X <- matrix(rlaplace(N * P), ncol=P)
    Y <- matrix(rmvnorm(N, mean=rep(0, P)), ncol=P)
    opt <- optim(theta.init, fn=J, method="CG", X=X, Y=Y,
		 control=list(maxit=1e3))
    print(c(opt$par, opt$value, mean(theta.p + opt$par)^2))
}


