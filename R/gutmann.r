library(Rcpp)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("gutmann.cpp")

N <- 150
A <- matrix(rnorm(16, 10, 10), nrow=4)
X <- matrix(rlaplace(4 * 150), nrow=N)
X <- X %*% A
Y <- matrix(rnorm(4*150), nrow=N)

opt <- optim(theta.init, J, method="CG", X=X, Y=Y)
theta.p <- c(t(A), log(abs(det(solve(A)))/4))

log10(mean((theta.p - opt$par)^2))

