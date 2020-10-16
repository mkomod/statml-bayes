#
# Replication of Gutmann & Hyvarinen Simulation
#
# ------------------------------------------------------------------------------
# Title: Noise-contrastive estimation
# Author: Gutmann, Hyvarinen
# Year: 2010
# ------------------------------------------------------------------------------
#
library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp("./rlaplace.cpp")
Rcpp::sourceCpp("./gutman.cpp")


# Settings ---------------------------------------------------------------------
set.seed(123)
N <- 10^2.4                            # Number of data points
P <- 2                                 # Dimensionality of data
Sig <- 1
mse <- function(theta, theta.p) mean((theta - theta.p) ^ 2)


# Data -------------------------------------------------------------------------
A <- matrix(rnorm(P^2, -10, 10), nrow=P, byrow=T)
X <- t(sapply(1:N, function(tmp) {
    s <- rlaplace(P)
    return( A %*% s)
}))
Y <- mvtnorm::rmvnorm(N, mean = rep(0, P), sigma=diag(Sig, P))


# Estimator --------------------------------------------------------------------
p.n <- function(x) mvtnorm::dmvnorm(x, rep(0, P), diag(Sig, P))

p.m <- function(x, B, cs) {
    ln.pm.0 <- -sqrt(2) * sum(abs(B %*% x))
    return(exp(ln.pm.0 + cs))
}

h <- function(x, B, cs) {
    p.m <- p.m(x, B, cs)
    p.n <- p.n(x)
    return( p.m / (p.m + p.n) )
}

J <- function(theta) {
    B <- solve(matrix(theta[1:(P^2)], nrow=P, byrow=T))
    cs <- theta[P^2 + 1]
    T.x <- apply(X, 1, function(x) { log(h(x, B, cs)) })
    T.y <- apply(Y, 1, function(y) { 1 - log(h(y, B, cs)) })
    return(1 / (2*N) * sum(T.x + T.y))
}


# Optimisation ----------------------------------------------------------------
J(theta.p)
theta.init <- rnorm(P^2 + 1)
opt <- optim(theta.init, fn=J, method="CG", # Conjugate Gradient
    control=list(fnscale=-1, reltol=1e-7))


# Results ---------------------------------------------------------------------
theta.p  <- c(t(A), log(abs(det(A))/4))
theta.hat <- opt$par
log10(mse(theta.hat, theta.p))

