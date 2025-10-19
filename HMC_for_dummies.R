
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


### Hamiltonian dynamics illustration

pdf("hamiltnoian_dynamics1.pdf", width = 12, height = 5)

par(mfrow = c(1,2))
r <- 3
a <- 1
t <- seq(0, 10, length = 5e2)
q <- r * cos(a + t)
p <- -r * sin(a + t)


plot(q, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3), asp = 1, lty = 2)
lines(q[1:101], p[1:101])
points(q[1], p[1], col = "blue", pch = 16)
text(q[1]+.8, p[1] - .05, "t = 0")

points(q[101], p[101], col = "blue", pch = 16)
text(q[101]+.8, p[101] - .05, "t = 2")
qt <- seq(-3, 3, length = 1e3)

plot(q, -dnorm(q, log = TRUE), type = 'l', ylab = "U(q)")
lines(q[1:101], -dnorm(q[1:101], log = TRUE)+ .06, col = "blue", lwd = 2)
points(x = q[1], y = -dnorm(q[1], log = TRUE), pch = 16, col = "blue")
text(q[1]+.8, -dnorm(q[1], log = TRUE), "t = 0")
points(x = q[101], y = -dnorm(q[101], log = TRUE) + .06, pch = 16, col = "blue")

text(q[101]+.8, -dnorm(q[101], log = TRUE), "t = 2")

dev.off()


pdf("hamiltnoian_dynamics2.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
q_new <- q
p_new <- p
p_new[1] <- 1

r_new <- -sqrt(q_new[1]^2 + p_new[1]^2)
a_new <-  acos(q_new[1]/r_new)

t <- seq(0, 10, length = 5e2)
q_new <- r_new * cos(a_new + t)
p_new <- -r_new * sin(a_new + t)


plot(q, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3), asp = 1, lty = 2)
lines(q[1:101], p[1:101])
points(q[1], p[1], col = "blue", pch = 16)

lines(q_new, p_new, col = "purple", lty = 2)
lines(q_new[1:101], p_new[1:101], col = "purple", lty = 1)
points(q_new[1], p_new[1], col = "orange", pch = 16)
text(q_new[1]+.8, p_new[1] - .05, "t = 0")
points(q[101], p[101], col = "blue", pch = 16)

points(q_new[101], p_new[101], col = "orange", pch = 16)
text(q_new[101]+.8, p_new[101] - .05, "t = 2")


qt <- seq(-3, 3, length = 1e3)

plot(q, -dnorm(q, log = TRUE), type = 'l', ylab = "U(q)")
lines(q_new[1:101], -dnorm(q_new[1:101], log = TRUE)+ .08, col = "purple", lwd = 2)
points(x = q_new[1], y = -dnorm(q_new[1], log = TRUE), pch = 16, col = "orange")
text(q_new[1]+.8, -dnorm(q_new[1], log = TRUE), "t = 0")
points(x = q_new[101], y = -dnorm(q_new[101], log = TRUE) + .1, pch = 16, col = "orange")
text(q_new[101]+.6, -dnorm(q_new[101], log = TRUE), "t = 2")

dev.off()


##### Gaussian HMC

normalHMC <- function(s = 1, n = 1e4)
{
  qt <- numeric(length = n)
  qt[1] <- 0 # starting here
  # don't need a starting value of p
  
  for(k in 2:n)
  {
    # new momentum
    p <- rnorm(1)
    
    #initial conditions for new H
    r2 <- qt[k-1]^2 + p^2
    # choosing +- with probability 1/2
    r <- sample(c(sqrt(r2), -sqrt(r2)), size = 1)  
    a <- acos(qt[k-1]/r)
    
    # simulating Hamiltonian forward s time units
    qt[k] <- r*cos(a + s)
    p <- -r * sin(a + s)
  }
  return(qt)
}

chain1 <- normalHMC(s = .1)
chain2 <- normalHMC(s = 1)
chain3 <- normalHMC(s = 5)

pdf("normal_hmc.pdf", width = 12, height = 5)

par(mfrow = c(1,2))
foo <- seq(-3, 3, length = 1e3)
plot(foo, dnorm(foo), type = 'l', ylab = "Density")
lines(density(chain3), col = "orange")
lines(density(chain2), col = "red")
lines(density(chain1), col = "blue")
legend("topright", 
       col = c("black", "orange", "red", "blue"), 
       legend = c("Truth", "s = 5", "s = 1", "s = .1"), lty = 1)

plot.ts(chain3, col = "orange", ylab = "Trace Plot", ylim = c(-4, 6))
lines(chain2, col = "red")
lines(chain1, col = "blue")
legend("topright", 
       col = c("orange", "red", "blue"), 
       legend = c("s = 5", "s = 1", "s = .1"), lty = 1)

dev.off()


#### Leaf-frog discretization

pdf("lfd1.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
foo <- seq(-3, 3, length = 1e3)
plot(foo, dnorm(foo), type = 'l', ylab = "Density")
lines(density(chain3), col = "orange")
lines(density(chain2), col = "red")
lines(density(chain1), col = "blue")
legend("topright", 
       col = c("black", "orange", "red", "blue"), 
       legend = c("Truth", "s = 5", "s = 1", "s = .1"), lty = 1)

plot.ts(chain3, col = "orange", ylab = "Trace Plot", ylim = c(-4, 6))
lines(chain2, col = "red")
lines(chain1, col = "blue")
legend("topright", 
       col = c("orange", "red", "blue"), 
       legend = c("s = 5", "s = 1", "s = .1"), 
       lty = 1)

dev.off()


pdf("lfd5.pdf", width = 4, height = 4)
r <- 3
a <- 1
t <- seq(0, 10, length = 5e2)
q <- r * cos(a + t)
p <- -r * sin(a + t)

eps <- .5
t <- seq(0,20*eps, by = eps)
e_p <- numeric(length = length(t))
e_q <- numeric(length = length(t))
e_p[1] <- p[1]
e_q[1] <- q[1]
for(i in 2:length(t))
{
  e_p[i] <- e_p[i-1] - eps*e_q[i-1]
  e_q[i] <- e_q[i-1] + eps*e_p[i]
}

plot(q, p, type = 'l', ylim = c(-3.5,3.5), 
     xlim = c(-3.5,3.5), asp = 1, lty = 1)
lines(e_q, e_p, lty = 2, col = "purple", type = "b", pch = 16)

eps <- .2
t <- seq(0,20*eps, by = eps)
e_p <- numeric(length = length(t))
e_q <- numeric(length = length(t))
e_p[1] <- p[1]
e_q[1] <- q[1]
for(i in 2:length(t))
{
  e_p[i] <- e_p[i-1] - eps*e_q[i-1]
  e_q[i] <- e_q[i-1] + eps*e_p[i]
}
lines(e_q, e_p, lty = 2, col = "blue", type = "b", pch = 16)

eps <- .1
t <- seq(0,20*eps, by = eps)
e_p <- numeric(length = length(t))
e_q <- numeric(length = length(t))
e_p[1] <- p[1]
e_q[1] <- q[1]
for(i in 2:length(t))
{
  e_p[i] <- e_p[i-1] - eps*e_q[i-1]
  e_q[i] <- e_q[i-1] + eps*e_p[i]
}
lines(e_q, e_p, lty = 2, col = "pink", type = "b", pch = 16)
legend("center", 
       legend = c("eps = .5", "eps = .2", "eps = .1"), 
       fill = c("purple", "blue", "pink"), cex = .80, bty = "n")
dev.off()



pdf("lfd4.pdf", width = 4, height = 4)
r <- 3
a <- 1
t <- seq(0, 10, length = 5e2)
q <- r * cos(a + t)
p <- -r * sin(a + t)
L <- 20

eps <- .3
e_p <- numeric(length = length(L))
e_q <- numeric(length = length(L))
e_q[1] <- q[1]
e_p[1] <- p[1] - eps/2 * e_q[1]
for(i in 2:L)
{
  e_q[i] <- e_q[i-1] + eps*e_p[i-1]
  if(i != L) e_p[i] <- e_p[i-1] - eps*e_q[i]
}
e_p[L] <- e_p[L-1] - eps*e_q[L]/2
plot(q, p, type = 'l', ylim = c(-3.5,3.5), xlim = c(-3.5,3.5), asp = 1, lty = 1)
lines(e_q, e_p, lty = 2, col = "purple", type = "b", pch = 16)

eps <- .2
e_p <- numeric(length = length(L))
e_q <- numeric(length = length(L))
e_q[1] <- q[1]
e_p[1] <- p[1] - eps/2 * e_q[1]
for(i in 2:L)
{
  e_q[i] <- e_q[i-1] + eps*e_p[i-1]
  if(i != L) e_p[i] <- e_p[i-1] - eps*e_q[i]
}
e_p[L] <- e_p[L-1] - eps*e_q[L]/2
lines(e_q, e_p, lty = 2, col = "blue", type = "b", pch = 16)

eps <- .1
e_p <- numeric(length = length(L))
e_q <- numeric(length = length(L))
e_q[1] <- q[1]
e_p[1] <- p[1] - eps/2 * e_q[1]
for(i in 2:L)
{
  e_q[i] <- e_q[i-1] + eps*e_p[i-1]
  if(i != L) e_p[i] <- e_p[i-1] - eps*e_q[i]
}
e_p[L] <- e_p[L-1] - eps*e_q[L]/2
lines(e_q, e_p, lty = 2, col = "pink", type = "b", pch = 16)
legend("center", 
       legend = c("eps = .3", "eps = .2", "eps = .1"), 
       fill = c("purple", "blue", "pink"), 
       cex = .80, bty = "n")

dev.off()


normalLF_HMC <- function(L = 10, eps = .1, n = 1e4)
{
  qt <- numeric(length = n)
  qt[1] <- 0 # starting here
  accept <- 1
  
  # vectorizing this outside to save time
  momentums <- rnorm(n) 
  
  for(k in 2:n)
  {
    # new momentum
    p <- momentums[k]
    q_prop <- qt[k-1]
    
    # one half-Euler step
    p_prop <- p - eps/2 * q_prop  
    for(l in 1:L) # let the frog leap!
    {
      q_prop <- q_prop + eps * p_prop
      if(l != L) p_prop <- p_prop - eps*q_prop
    }
    # one last half-Euler step
    p_prop <- p_prop - eps/2 * q_prop
    
    # Accept-reject
    log.ratio <-  (-q_prop^2 - p_prop^2 + qt[k-1]^2 + p^2)/2
    if(log(runif(1)) < log.ratio)
    {
      qt[k] <- q_prop
      accept <- accept + 1
    } else{
      qt[k] <- qt[k-1]
    }
  }
  print(paste("Acceptance = ", accept/n))
  return(qt)
}

#Keeping L * epsilon = 1
chain1 <- normalLF_HMC(L = 10, eps = .1)
chain2 <- normalLF_HMC(L = 1, eps = 1)
chain3 <- normalLF_HMC(L = 100, eps = .01) # more time consuming

pdf("lfd2.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
foo <- seq(-3, 3, length = 1e3)
plot(foo, dnorm(foo), type = 'l', 
     ylab = "Density",
     xlab = "q")
lines(density(chain3), col = "orange")
lines(density(chain2), col = "red")
lines(density(chain1), col = "blue")
legend("topright", 
       col = c("black", "orange", "red", "blue"), 
       legend = c("Truth", "eps = .01", "eps = 1", "eps = .1"), 
       lty = 1)

plot.ts(chain3, col = "orange", ylab = "Trace Plot", ylim = c(-4, 6))
lines(chain2, col = "red")
lines(chain1, col = "blue")
legend("topright", 
       col = c("orange", "red", "blue"), 
       legend = c("eps = .01", "eps = 1", "eps = .1"), 
       lty = 1)
dev.off()


#Keeping L * epsilon = 10
chain1 <- normalLF_HMC(L = 100, eps = .1) # more time consuming
chain2 <- normalLF_HMC(L = 10, eps = 1)
chain3 <- normalLF_HMC(L = 1, eps = 10) #cheap but inaccurate

pdf("lfd3.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
foo <- seq(-3, 3, length = 1e3)
plot(foo, dnorm(foo), type = 'l',
     ylab = "Density",
     xlab = "q")
lines(density(chain2), col = "red")
lines(density(chain1), col = "blue")
lines(density(chain3), col = "orange")
legend("topright", 
       col = c("black", "orange", "red", "blue"), 
       legend = c("Truth", "eps = 10", "eps = 1", "eps = .1"), 
       lty = 1)

plot.ts(chain2, col = "red", ylab = "Trace Plot", ylim = c(-4, 6))
lines(chain1, col = "blue")
lines(chain3, col = "orange")

legend("topright", col = c("orange", "red", "blue"), 
       legend = c("eps = 10", "eps = 1", "eps = .1"), 
       lty = 1)

dev.off()
