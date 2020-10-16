library(Rcpp)
library(parallel)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("gutmann.cpp")

CORES <- 20
NUM_DATASETS <- 600
NUM_START_VALS <- 10

set.seed(2020)

generate_data <- function(N) {
    A <- matrix(c(rnorm(4, 0, 0.1), rnorm(4, 3, 0.1),
		  rnorm(4, -2, 0.1), rnorm(4, 1, 0.1)), nrow=4)
    X <- matrix(rlaplace(4 * N), nrow=N) %*% A
    Y <- matrix(rnorm(4 *  N), nrow=N)
    return(list(X=X, Y=Y, A=A))
}

vals <- mclapply(seq(2, 4, length.out=CORES), function (samp.exp) {
    res <- matrix(0, nrow=NUM_DATASETS, ncol=NUM_START_VALS)

    for (i in 1:NUM_DATASETS) {
	# Generate the data for this run
	l <- generate_data(10^samp.exp)
	X <- l$X; Y <- l$Y; A <- l$A
	for (j in 1:NUM_START_VALS) {
	    theta.init <- c(rnorm(16, 1, sd=2.2), runif(1, -10, -1))
	    tryCatch({
		opt <- optim(theta.init, J, method="CG", X=X, Y=Y)
		theta.p <- c(t(A), log(abs(det(solve(A)))/4))
		res[i, j] <- mean((theta.p - opt$par)^2)
		cat(".")
	    }, error = function(e) {
		# On error set value to NA
		print(e)
		res[i, j] <- NA
		cat("x")
	    })
	}
    }
    res[res == 0] <- NA
    return(res)
}, mc.cores=CORES)

save(vals, file="vals.RData")

