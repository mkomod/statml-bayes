

NS <- c(100, 193, 373, 720, 1389, 2683, 3162, 5179, 10000)

N <- NS[1]
res <- data.frame(lwr=numeric(), med=numeric(), upr=numeric())

for (N in NS) {
load(paste0("./m_output/r_", N, ".RData"))
r <- get(paste0("r_", N))
vals <- by(r, r$run, function(x) min(x$mse))
rs <- quantile(vals, c(0.3, 0.5, 0.7))
res <- rbind(res, list(lwr=rs[1], med=rs[2], upr=rs[3]))
}

plot(log10(res$med), ylim=c(-2.2, -0.5), type="b", pch=16)
x <- 1:9
polygon(c(x, rev(x)), log10(c(res$upr, rev(res$lwr))), border=NA, col=rgb(0, 1, 0, 0.3))
