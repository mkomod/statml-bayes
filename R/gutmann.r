library(Rcpp)
library(parallel)
library(mvtnorm)

Rcpp::sourceCpp("rlaplace.cpp")
Rcpp::sourceCpp("gutmann.cpp")

CORES <- 1
NUM_DATASETS <- 200
NUM_START_VALS <- 5

set.seed(2020)


generate_data <- function(N) {
    E <- TRUE
    while (E == TRUE) {
	tryCatch({
	    A <- matrix(round(c(runif(4, 4, 6), runif(4, 6, 8),
				runif(4, 4, 6), runif(4, 3, 5)), 1), nrow=4)
	    X <- matrix(rlaplace(4 * N), nrow=N) %*% t(A)
	    S <- A %*% t(A)
	    Y <- matrix(rmvnorm(N, rep(0, 4), S), nrow=N, byrow=T)
	    E <- FALSE
	}, error = function(e) { E <- TRUE })
    }	

    return(list(X=X, Y=Y, A=A, S=S))
}

generate_intial_theta <- function() {
    return(c(
	     round(runif(16, 3, 8), 2), 
	     runif(1, -5, -2.5)
    ))
}


vals <- mclapply(seq(3, 3, length.out=CORES), function (samp.exp) {
    samp.exp <- round(samp.exp)        # round sample size to int
    res <- matrix(0, nrow=NUM_DATASETS, ncol=NUM_START_VALS)

    for (i in 1:NUM_DATASETS) {
	# Generate the data for this run
	l <- generate_data(10^samp.exp)
	X <- l$X; Y <- l$Y; A <- l$A; S <- l$S   # Unpack list

	# Needed for normal density
	S_inv <- solve(S);
	lNC <- log(sqrt((2 * pi)^4 * det(S) ))
	
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



