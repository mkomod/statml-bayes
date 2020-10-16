library(Rcpp)
library(parallel)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("gutmann.cpp")

CORES <- 4

vals <- mclapply(seq(2, 4, length.out=CORES), function (samp.exp) {
    val <- c()
    for (i in 1:700) {
	l <- generate_data(10^samp.exp)
	X <- l$X; Y <- l$Y
	for (j in 1:5) {
	    theta.init <- rnorm(17, sd=3)
	    tryCatch({
		opt <- optim(theta.init, J, method="CG", X=X, Y=Y)
		theta.p <- c(t(A), log(abs(det(solve(A)))/4))
		val[i*(j-1) + i] <- mean((theta.p - opt$par)^2)
	    }, error = function(e) {
		val[i*(j-1) + i] <- NA
	    })
	}
    }
    return(val)
}, mc.cores=CORES)

generate_data <- function(N) {
    A <- matrix(c(rnorm(4, 0, 0.1), rnorm(4, 3, 0.1), 
		  rnorm(4, 0, 0.1), rnorm(4, -3, 0.1)) , nrow=4)
    X <- matrix(rlaplace(4 * N), nrow=N) %*% A
    Y <- matrix(rnorm(4 *  N), nrow=N)
    return(list(X=X, Y=Y))
}

