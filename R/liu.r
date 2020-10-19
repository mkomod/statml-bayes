
set.seed(123)
discriminance_sampling <- function(z, X) {
    p0 <- dnorm(X, 0, s0)
    1/2 - 1 / (2 * N) * sum(p0 / (dnorm(X, 0, 1) / z + p0))
}

importance_sampling <- function(X) {
    mean(dnorm(X, 0, 1) / dnorm(X, 0, s0) )
}

reverse_importance_sampling <- function(X) {
    1 / mean(dnorm(X, 0, s0) / dnorm(X, 0, 1)) 
}

N <- 1e4; NN <- 50
S0 <- c(seq(from=1/4, to=1, length.out=15), 
	seq(from=1.1, to=4, length.out=15))
m.ds <- m.is <- m.ris <- matrix(nrow=NN, ncol=length(S0))
for (i in 1:NN) {
    val.ds <- val.is <- val.ris <- c()
    for (s0 in S0) {
	X0 <- rnorm(N, 0, 1)
	X1 <- rnorm(N, 0, s0)
	X <- c(X0, X1)
	val.ds <- c(val.ds, uniroot(discriminance_sampling, c(0, 100), X=X)$root)
	val.is <- c(val.is, importance_sampling(X1))
	val.ris <- c(val.ris, reverse_importance_sampling(X0))
    }
    m.ds[i, ] <- val.ds
    m.is[i, ] <- val.is
    m.ris[i, ] <- val.ris
}
m.m.ds <- apply(m.ds, 2, mean)
m.m.is <- apply(m.is, 2, mean)
m.m.ris <- apply(m.ris, 2, mean)

pdf("liu.pdf", 9, 6)
par(mar=c(5, 5, 1, 1))
plot(S0, log(m.m.ris), type="b", col="green", pch=16,
     ylim=c(-0.4, 0.4), log="x", 
     ylab="logZ", xlab=expression(sigma[0]))
points(S0, log(m.m.is), type="b", col="red", pch=17)
points(S0, log(m.m.ds), type="b", col="blue", pch=18)
legend("bottomright", 
       legend=c("RIS", "IS", "DS"), 
       lty=c(1, 1, 1), pch=c(16, 17, 18), col=c("green","red", "blue"))
dev.off()
