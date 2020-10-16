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

Rcpp::sourceCpp("./rlaplace.cpp")


# Settings ---------------------------------------------------------------------
N <- 10^2.4                            # Number of data points
P <- 4                                 # Dimensionality of data
Sig <- 1


# Data -------------------------------------------------------------------------
A <- diag(Sig, P)                        # Covariance Matrix
X <- t(sapply(1:N, function(tmp) {
    s <- rlaplace(P)
    return( A %*% s)
}))
Y <- mvtnorm::rmvnorm(N, mean = rep(0, P), sigma=A)


# Estimator --------------------------------------------------------------------
p.n <- function(x) mvtnorm::dmvnorm(x, rep(0, P), A)

p.m <- function(x, B, cs) {
    ln.pm.0 <- -sqrt(2) * sum(abs(t(x) %*% t(B)))
    return(exp(ln.pm.0 + cs))
}

h <- function(x, B, c) {
    p <- p.m(x, B, c)
    q <- p.n(x)
    return( p / (p + q) )
}

J <- function(theta) {
    if (theta[1] < 0) {
	T.y <- T.x <- -Inf
    } else {
	B <- solve(diag(theta[1], P))
	cs <- theta[2]
	T.x <- apply(X, 1, function(x) { log(h(x, B, cs)) })
	T.y <- apply(Y, 1, function(y) { 1 - log(h(y, B, cs)) })
    }
    return(1 / (2*N) * sum(T.x + T.y))
}


# Optimisation ----------------------------------------------------------------
set.seed(201014)
theta.init <- c(abs(rnorm(1, sd=3)), rnorm(1, sd=3)) # Initial theta values
opt <- optim(theta.init, fn=J,
	     method="CG",                     # Conjugate Gradient
	     control=list(fnscale=-1,         # fnscale=1 for max (min by def)
			  reltol=1e-4,
			  type=3)          
)


# Results ---------------------------------------------------------------------
mse <- function(theta, theta.p) mean((theta - theta.p) ^ 2)
theta.p  <- c(rep(Sig, P), log(abs(det(A))/4))
opt.scl <- opt$par / opt$par[1]
theta.hat <- c(rep(opt.scl[1], P), opt.scl[2])

mse(theta.hat, theta.p) 
log10(mse(theta.hat, theta.p))

