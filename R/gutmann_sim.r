# Replication of Gutmann & Hyvarinen Simulation
#
# ------------------------------------------------------------------------------
# Title: Noise-contrastive estimation
# Author: Gutmann, Hyvarinen
# Year: 2010
# ------------------------------------------------------------------------------
#
library(Rcpp)
library(mvtnorm)
library(numDeriv)

Rcpp::sourceCpp("./rlaplace.cpp")


# Settings ---------------------------------------------------------------------
N <- 10^2.4                            # Number of data points
P <- 4                                 # Dimensionality of data


# Data -------------------------------------------------------------------------
A <- diag(1, P)                        # Covariance Matrix
X <- t(sapply(1:N, function(tmp) {
    s <- rlaplace(P)
    A %*% s
}))
Y <- mvtnorm::rmvnorm(N, mean = rep(0, P), sigma=A)


# Estimator --------------------------------------------------------------------
p.n <- function(x) mvtnorm::dmvnorm(x, rep(0, P), A)

p.m <- function(x, theta) {
    A <- matrix(theta[1:(P^2)], nrow = P, byrow=TRUE)
    B <- solve(A)
    c <- theta[P^2 + 1]
    ln.pm.0 <- -sum(apply(B, 1, function(b) sqrt(2) * abs(b %*% x)))
    return(exp(ln.pm.0 + c))
}

h <- function(x, theta) {
    return(p.m(x, theta) / (p.m(x, theta) + p.n(x)))
}

J <- function(theta) {
    T.x <- apply(X, 1, function(x) { log(h(x, theta)) })
    T.y <- apply(Y, 1, function(y) { 1 -log(h(y, theta)) })
    return(1 / (2*N) * sum(T.x + T.y))
}

J.Delta <- function(theta) {           # Gradient of J
    return(numDeriv::grad(J, theta))
}


# Optimisation ----------------------------------------------------------------
set.seed(201014)
theta.init = c(runif(P^2), 0)          # Initial theta values
opt <- optim(theta.init, fn=J, gr=J.Delta, 
     method="CG",                      # Conjugate Gradient
     control=list(fnscale=-1)          # fnscale=1 for max (min by def)
)


# Results ---------------------------------------------------------------------
mse <- function(theta, theta.p) mean((theta - theta.p) ^ 2)
theta.p  <- c(as.vector(A), log(abs(det(solve(A))) / 4))
theta.hat <- opt$par

log10(mse(theta.hat, theta.p))

