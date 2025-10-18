
mydir <- this.path::here()
setwd(mydir)


### HMC for dummies

pdf("gaussian_potential.pdf", width = 8, height = 4)
q <- seq(-3, 3, length = 1e3)
par(mfrow = c(1,2))
plot(q, dnorm(q), type = 'l', ylab = "Density of N(0,1)")
points(x = -2, y = 0, pch = 16, col = "blue")
plot(q, -dnorm(q, log = TRUE), type = 'l', ylab = "Negative log-density of N(0,1)")
points(x = -2, y = -dnorm(-2, log = TRUE) + .1, pch = 16, col = "blue")

dev.off()
