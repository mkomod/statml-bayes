library(Rcpp)
library(parallel)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("gutmann.cpp")

CORES <- 6
NUM_DATASETS <- 200
NUM_START_VALS <- 5

set.seed(2020)


generate_data <- function(N) {
    A <- matrix(round(runif(16, 0, 10), 2), nrow=4)
    X <- matrix(rlaplace(4 * N), nrow=N) %*% A
    Y <- matrix(rnorm(4 *  N), nrow=N)
    return(list(X=X, Y=Y, A=A))
}

generate_intial_theta <- function() {
    return(c(
	     round(runif(16, 3, 8), 1), 
	     runif(1, -8, -1)
    ))
}


vals <- mclapply(seq(2, 3.5, length.out=CORES), function (samp.exp) {
    samp.exp <- round(samp.exp)        # round sample size to int
    res <- matrix(0, nrow=NUM_DATASETS, ncol=NUM_START_VALS)

    for (i in 1:NUM_DATASETS) {
	# Generate the data for this run
	E <- TRUE
	while (E == TRUE) {
	    tryCatch({
		l <- generate_data(10^samp.exp)
		X <- l$X; Y <- l$Y; A <- l$A          # Unpack list

		# Needed for normal density
		S <- A %*% t(A)
		S_inv <- solve(S);
		lNC <- log(sqrt((2 * pi)^4 * det(S) ))
		E <- FALSE
	    }, error = function(e) E <- TRUE)
	}	

	# Theta values used to generate the data and Normalising constant
	theta.p <- c(t(A), log(abs(det(solve(A)))/4))

	# Solve for different starting values and save MSE
	for (j in 1:NUM_START_VALS) {
	    tryCatch({
		theta.init <- generate_intial_theta()
		opt <- optim(theta.init, J, 
		    method="CG", 
		    X=X, 
		    Y=Y, 
		    S_inv=S_inv, 
		    log_NormalisingConst=lNC
		)
		res[i, j] <- mean((theta.p - opt$par)^2)
		cat(".")
	    }, error = function(e) {
		# print(e)
		res[i, j] <- NA	                     # Sets value to 0 not NA
		cat("x")
	    })
	}

    }
    res[res == 0] <- NA                # Set 0s to NA
    return(res)

}, mc.cores=CORES)


save(vals, file="vals.RData")          # Save data for later use

