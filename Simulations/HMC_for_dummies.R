rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

pdf("gaussian_potential.pdf", width = 12, height = 5)

x <- seq(-3, 3, length = 5e2)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

# Left plot: Density
plot(x, dnorm(x), type = 'l',
     ylab = expression(Density~pi(x)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2)   # enlarge tick labels
points(x = -2, y = 0, pch = 16, col = "blue")

# Right plot: Potential
plot(x, -dnorm(x, log = TRUE), type = 'l',
     ylab = expression(Potential~U(x)),
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,
     cex.axis = 1.2)
points(x = -2, y = -dnorm(-2, log = TRUE) + .1, pch = 16, col = "blue")

dev.off()

### Hamiltonian dynamics illustration

pdf("hamiltnoian_dynamics1.pdf", width = 12, height = 5)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

r <- 3
a <- 1
t <- seq(0, 10, length = 5e2)
x <- r * cos(a + t)
p <- -r * sin(a + t)


plot(x, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3),
     asp = 1, lty = 2,  ylab = expression(Momentum~(p)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2)
lines(x[1:101], p[1:101])
points(x[1], p[1], col = "blue", pch = 16)
text(x[1]+.8, p[1] - .05, "t = 0")

points(x[101], p[101], col = "blue", pch = 16)
text(x[101]+.8, p[101] - .05, "t = 2")
qt <- seq(-3, 3, length = 1e3)

plot(x, -dnorm(x, log = TRUE), type = 'l',  
     ylab = expression(Potential~U(x)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2)
lines(x[1:101], -dnorm(x[1:101], log = TRUE)+ .06, col = "blue", lwd = 2)
points(x = x[1], y = -dnorm(x[1], log = TRUE), pch = 16, col = "blue")
text(x[1]+.8, -dnorm(x[1], log = TRUE), "t = 0")
points(x = x[101], y = -dnorm(x[101], log = TRUE) + .06, pch = 16, col = "blue")

text(x[101]+.8, -dnorm(x[101], log = TRUE), "t = 2")

dev.off()


pdf("hamiltnoian_dynamics2.pdf", width = 12, height = 5)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

x_new <- x
p_new <- p
p_new[1] <- 1

r_new <-  sqrt(x_new[1]^2 + p_new[1]^2)
a_new <-  acos(x_new[1]/r_new)

t <- seq(0, 10, length = 5e2)
x_new <- r_new * cos(a_new + t)
p_new <- -r_new * sin(a_new + t)


plot(x, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3), 
     asp = 1, lty = 2, ylab = expression(Momentum~(p)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2
     )
lines(x[1:101], p[1:101])
points(x[1], p[1], col = "blue", pch = 16)

lines(x_new, p_new, col = "purple", lty = 2)
lines(x_new[1:101], p_new[1:101], col = "purple", lty = 1)
points(x_new[1], p_new[1], col = "orange", pch = 16)
text(x_new[1]+.8, p_new[1] - .05, "t = 0")
points(x[101], p[101], col = "blue", pch = 16)

points(x_new[101], p_new[101], col = "orange", pch = 16)
text(x_new[101]+.8, p_new[101] - .05, "t = 2")


qt <- seq(-3, 3, length = 1e3)

plot(x, -dnorm(x, log = TRUE), type = 'l', 
     ylab = expression(Potential~(U(x))),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2)
lines(x_new[1:101], -dnorm(x_new[1:101], log = TRUE)+ .08, col = "purple", lwd = 2)
points(x = x_new[1], y = -dnorm(x_new[1], log = TRUE), pch = 16, col = "orange")
text(x_new[1]+.8, -dnorm(x_new[1], log = TRUE), "t = 0")
points(x = x_new[101], y = -dnorm(x_new[101], log = TRUE) + .1, pch = 16, col = "orange")
text(x_new[101]+.6, -dnorm(x_new[101], log = TRUE), "t = 2")

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

pdf("exact_normal1.pdf", width = 16, height = 5)

# layout and base margins
par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))

# desired sizes (tweak here)
cex_lab  <- 1.5   # axis label size
cex_axis <- 1.2   # tick label size
leg_cex  <- 2

x <- seq(-3, 3, length = 1e3)

## (1) Density plot: explicitly supply cex.* here so it matches
plot(x, dnorm(x), type = "l", lwd = 2,
     ylab = expression(Density~pi(x)), xlab = expression(Position~(x)),
     cex.lab = cex_lab, cex.axis = cex_axis)
lines(density(chain3), col = "orange", lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain1), col = "blue",   lwd = 2)
legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth", "s = 5", "s = 1", "s = .1"),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot: suppress axes/labels and add them manually with cex
plot.ts(chain3, col = "orange", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE,  cex.lab = cex_lab)   # suppress default axes and labels
lines(chain2, col = "red",  lwd = 2)
lines(chain1, col = "blue", lwd = 2)
# add axes with explicit cex.axis
axis(1, cex.axis = cex_axis)   # x-axis
axis(2, cex.axis = cex_axis)   # y-axis
box()
# add axis labels with explicit cex
mtext(expression(Iterations), side = 1, line = 2.8, cex = cex_lab)
mtext(expression("Trace plot"), side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("orange", "red", "blue"),
       legend = c("s = 5", "s = 1", "s = .1"),
       lty = 1, cex = 2, bty = "n")

## (3) ACF as line plots: compute acf objects then plot manually
lagmax <- 50
acf1 <- acf(chain1, plot = FALSE, lag.max = lagmax)
acf2 <- acf(chain2, plot = FALSE, lag.max = lagmax)
acf3 <- acf(chain3, plot = FALSE, lag.max = lagmax)

# ensure lags are numeric vector and same length
lags <- as.numeric(acf1$lag)

# plot first series (no axes)
plot(lags, as.numeric(acf1$acf), type = "l", lwd = 2, col = "blue",
     ylim = c(min(c(acf1$acf, acf2$acf, acf3$acf)),
              max(c(acf1$acf, acf2$acf, acf3$acf))),
     axes = FALSE, xlab = "", ylab = "", cex.lab = cex_lab)
lines(lags, as.numeric(acf2$acf), col = "red",    lwd = 2)
lines(lags, as.numeric(acf3$acf), col = "orange", lwd = 2)

# add axes with explicit sizes
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext("Lag", side = 1, line = 2.8, cex = cex_lab)
mtext("ACF", side = 2, line = 3.5, cex = cex_lab)

legend("topright",
       col = c("orange", "red", "blue"),
       legend = c("s = 5", "s = 1", "s = .1"),
       lty = 1, cex = 2, bty = "n")

dev.off()



#### Leaf-frog discretization with different s

pdf("lfd_with_diff_s.pdf", width = 12, height = 5)
par(mfrow = c(1,2))
x <- seq(-3, 3, length = 1e3)
plot(x, dnorm(x), type = 'l',  lwd = 2,
     ylab = expression(Density~pi(x)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2)
lines(density(chain3),  lwd = 2, col = "orange")
lines(density(chain2),  lwd = 2, col = "red")
lines(density(chain1),  lwd = 2, col = "blue")
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

## Euler integrator — one panel per epsilon
pdf("eular_traj.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0,0,0,0))

# Plotting settings
cex_lab  <- 1.5
cex_axis <- 1.2
cex_leg  <- 2
pt_cex   <- 1.5
lwd_line <- 2

# True Hamiltonian (analytic circular orbit)
r <- 3
a <- 1
t_cont <- seq(0, 2*pi, length.out = 1e3)
x_true <- r * cos(a + t_cont)
p_true <- -r * sin(a + t_cont)

eps_vals <- c(0.5, 0.2, 0.1)
cols <- c("darkgreen", "darkblue", "darkred")

for (k in seq_along(eps_vals)) {
  eps <- eps_vals[k]
  
  # Euler discrete trajectory
  t_seq <- seq(0, 20 * eps, by = eps)
  n <- length(t_seq)
  e_p <- numeric(n)
  e_q <- numeric(n)
  e_p[1] <- p_true[1]
  e_q[1] <- x_true[1]
  for (i in 2:n) {
    e_p[i] <- e_p[i-1] - eps * e_q[i-1]
    e_q[i] <- e_q[i-1] + eps * e_p[i-1]
  }
  
  # Panel plot: true contour + Euler trajectory
  plot(x_true, p_true, type = "l", lwd = lwd_line, col = "black",
       xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), asp = 1,
       xlab = expression(Position~(x)),
       ylab = expression(Momentum~(p)),
       cex.lab = cex_lab, cex.axis = cex_axis, las = 1)
  lines(e_q, e_p, type = "b", pch = 16, cex = pt_cex,
        col = cols[k], lty = 2, lwd = lwd_line)
  
  # Legend: true Hamiltonian + epsilon
  legend("topright",
         legend = c(expression(H_true), bquote(epsilon == .(eps))),
         col = c("black", cols[k]),
         lty = c(1,2),
         pch = c(NA, 16),
         lwd = c(lwd_line, lwd_line),
         bty = "n", cex = cex_leg)
}

dev.off()

## Modified Euler integrator — one panel per epsilon
pdf("modified_euler_traj.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0,0,0,0))

# Plotting settings
cex_lab  <- 1.5
cex_axis <- 1.2
cex_leg  <- 2
pt_cex   <- 1.5
lwd_line <- 2

# True Hamiltonian (analytic circular orbit)
r <- 3
a <- 1
t_cont <- seq(0, 2*pi, length.out = 1e3)
x_true <- r * cos(a + t_cont)
p_true <- -r * sin(a + t_cont)

eps_vals <- c(0.5, 0.2, 0.1)
cols <- c("darkgreen", "darkblue", "darkred")

for (k in seq_along(eps_vals)) {
  eps <- eps_vals[k]
  
  # Modified Euler discrete trajectory
  t_seq <- seq(0, 20 * eps, by = eps)
  n <- length(t_seq)
  m_p <- numeric(n)
  m_q <- numeric(n)
  m_p[1] <- p_true[1]
  m_q[1] <- x_true[1]
  
  for (i in 2:n) {
    # Modified Euler (symplectic Euler)
    m_p[i] <- m_p[i-1] - eps * m_q[i-1]
    m_q[i] <- m_q[i-1] + eps * m_p[i]   # <-- use new momentum!
  }
  
  # Panel plot: true contour + modified Euler trajectory
  plot(x_true, p_true, type = "l", lwd = lwd_line, col = "black",
       xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), asp = 1,
       xlab = expression(Position~(x)),
       ylab = expression(Momentum~(p)),
       cex.lab = cex_lab, cex.axis = cex_axis, las = 1)
  lines(m_q, m_p, type = "b", pch = 16, cex = pt_cex,
        col = cols[k], lty = 2, lwd = lwd_line)
  
  # Legend: true Hamiltonian + epsilon
  legend("topright",
         legend = c(expression(H_true), bquote(epsilon == .(eps))),
         col = c("black", cols[k]),
         lty = c(1,2),
         pch = c(NA, 16),
         lwd = c(lwd_line, lwd_line),
         bty = "n", cex = cex_leg)
}

dev.off()



## lfd6: Leapfrog integrator — one panel per epsilon
pdf("leapfrog_traj.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0,0,0,0))

# Plotting settings
cex_lab  <- 1.5
cex_axis <- 1.2
cex_leg  <- 2
pt_cex   <- 1.5
lwd_line <- 2

# True Hamiltonian (analytic circular orbit)
r <- 3
a <- 1
t_cont <- seq(0, 2*pi, length.out = 1000)
x_true <- r * cos(a + t_cont)
p_true <- -r * sin(a + t_cont)

eps_vals <- c(0.5, 0.2, 0.1)
cols <- c("darkgreen", "darkblue", "darkred")
L <- 20

for (k in seq_along(eps_vals)) {
  eps <- eps_vals[k]
  
  # Leapfrog integrator (position q, momentum p)
  e_q <- numeric(L)
  e_p <- numeric(L)
  e_q[1] <- x_true[1]
  e_p[1] <- p_true[1] - (eps/2) * e_q[1]
  if (L >= 2) {
    for (i in 2:L) {
      e_q[i] <- e_q[i-1] + eps * e_p[i-1]
      if (i != L) e_p[i] <- e_p[i-1] - eps * e_q[i]
    }
    e_p[L] <- e_p[L-1] - (eps/2) * e_q[L]
  }
  
  # Panel plot: true contour + leapfrog trajectory
  plot(x_true, p_true, type = "l", lwd = lwd_line, col = "black",
       xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), asp = 1,
       xlab = expression(Position~(x)),
       ylab = expression(Momentum~(p)),
       cex.lab = cex_lab, cex.axis = cex_axis, las = 1)
  lines(e_q, e_p, type = "b", pch = 16, cex = pt_cex,
        col = cols[k], lty = 2, lwd = lwd_line)
  
  # Legend: true Hamiltonian + epsilon
  legend("topright",
         legend = c(expression(H_true), bquote(epsilon == .(eps))),
         col = c("black", cols[k]),
         lty = c(1,2),
         pch = c(NA, 16),
         lwd = c(lwd_line, lwd_line),
         bty = "n", cex = cex_leg)
}

dev.off()


########## 

pdf("lfd_normal_diff_s.pdf", width = 4, height = 4)
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
lines(e_q, e_p, lty = 2, col = "darkgreen", type = "b", pch = 16)
legend("center", 
       legend = c("eps = .3", "eps = .2", "eps = .1"), 
       fill = c("darkgreen", "darkblue", "darkred"), 
       cex = .80, bty = "n")

dev.off()


########## Leapfrog for standard Normal ############


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

########## Leapfrog with s = 1 #########
pdf("lfd_with_s1.pdf", width = 16, height = 5)

# layout and base margins
par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))

# desired sizes (matching lfd1)
cex_lab  <- 1.5
cex_axis <- 1.2
leg_cex  <- 2

x <- seq(-3, 3, length = 1e3)

## (1) Density plot
plot(x, dnorm(x), type = "l", lwd = 2,
     ylab = expression(Density~pi(x)), xlab = expression(Position~(x)),
     cex.lab = 2.2, cex.axis = 1.8)
lines(density(chain3), col = "orange", lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain1), col = "blue",   lwd = 2)
legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth", 
                  expression(epsilon == 0.01),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot
plot.ts(chain3, col = "orange", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE)
lines(chain2, col = "red",  lwd = 2)
lines(chain1, col = "blue", lwd = 2)
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext(expression(Iterations), side = 1, line = 2.8, cex = cex_lab)
mtext("Trace plot", side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("orange", "red", "blue"),
       legend = c(expression(epsilon == 0.01),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2, bty = "n")

## (3) ACF as line plots
lagmax <- 50
acf1 <- acf(chain1, plot = FALSE, lag.max = lagmax)
acf2 <- acf(chain2, plot = FALSE, lag.max = lagmax)
acf3 <- acf(chain3, plot = FALSE, lag.max = lagmax)
lags <- as.numeric(acf1$lag)

plot(lags, as.numeric(acf1$acf), type = "l", lwd = 2, col = "blue",
     ylim = c(min(c(acf1$acf, acf2$acf, acf3$acf)),
              max(c(acf1$acf, acf2$acf, acf3$acf))),
     axes = FALSE, xlab = "", ylab = "")
lines(lags, as.numeric(acf2$acf), col = "red", lwd = 2)
lines(lags, as.numeric(acf3$acf), col = "orange", lwd = 2)
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext("Lag", side = 1, line = 2.8, cex = cex_lab)
mtext("ACF", side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("orange", "red", "blue"),
       legend = c(expression(epsilon == 0.01),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2, bty = "n")

dev.off()

# Keeping L * epsilon = 10
chain1 <- normalLF_HMC(L = 100, eps = .1)
chain2 <- normalLF_HMC(L = 10,  eps = 1)
chain3 <- normalLF_HMC(L = 1,   eps = 10)

########## Leapfrog with s = 10 #########
pdf("lfd_with_s10.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))

# desired sizes (matching lfd1)
cex_lab  <- 1.5
cex_axis <- 1.2
leg_cex  <- 2

x <- seq(-3, 3, length = 1e3)

## (1) Density plot
plot(x, dnorm(x), type = "l", lwd = 2,
     ylab = expression(Density~pi(x)), xlab = expression(Position~(x)),
     cex.lab = 2.2, cex.axis = 1.8)
lines(density(chain3), col = "orange", lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain1), col = "blue",   lwd = 2)
legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth",
                  expression(epsilon == 10),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot
plot.ts(chain3, col = "orange", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE)
lines(chain2, col = "red",  lwd = 2)
lines(chain1, col = "blue", lwd = 2)
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext(expression(Iterations), side = 1, line = 2.8, cex = cex_lab)
mtext("Trace plot", side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("orange", "red", "blue"),
       legend = c(expression(epsilon == 10),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2, bty = "n")

## (3) ACF as line plots
lagmax <- 50
acf1 <- acf(chain1, plot = FALSE, lag.max = lagmax)
acf2 <- acf(chain2, plot = FALSE, lag.max = lagmax)
acf3 <- acf(chain3, plot = FALSE, lag.max = lagmax)
lags <- as.numeric(acf1$lag)

plot(lags, as.numeric(acf1$acf), type = "l", lwd = 2, col = "blue",
     ylim = c(min(c(acf1$acf, acf2$acf, acf3$acf)),
              max(c(acf1$acf, acf2$acf, acf3$acf))),
     axes = FALSE, xlab = "", ylab = "")
lines(lags, as.numeric(acf2$acf), col = "red", lwd = 2)
lines(lags, as.numeric(acf3$acf), col = "orange", lwd = 2)
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext("Lag", side = 1, line = 2.8, cex = cex_lab)
mtext("ACF", side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("orange", "red", "blue"),
       legend = c(expression(epsilon == 10),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2, bty = "n")

dev.off()

#### Taking account for periodicity #######
eps <- 0.1

L_full    <- round(2 * pi / eps)        # ≈ 628
L_half    <- round(pi / eps)            # ≈ 314
L_quarter <- round(pi / (2 * eps))      # ≈ 157

cat("L values:\n")
cat("Full period L =", L_full, "\n")
cat("Half period L =", L_half, "\n")
cat("Quarter period L =", L_quarter, "\n")

########## Run HMC Chains ##########

chain_full    <- normalLF_HMC(L = L_full,    eps = eps)
chain_half    <- normalLF_HMC(L = L_half,    eps = eps)
chain_quarter <- normalLF_HMC(L = L_quarter, eps = eps)

########## Plot: Density, Trace, ACF ##########

pdf("lfd_period_based.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0,0,0,0))

cex_lab  <- 1.5
cex_axis <- 1.2
leg_cex  <- 2

x <- seq(-3, 3, length = 1000)

### (1) Density plot
plot(x, dnorm(x), type = "l", lwd = 2,
     ylab = expression(Density~pi(x)),
     xlab = expression(Position~(x)),
     cex.lab = 2.2, cex.axis = 1.8)

lines(density(chain_quarter), col = "orange", lwd = 2)
lines(density(chain_half),    col = "red",    lwd = 2)
lines(density(chain_full),    col = "blue",   lwd = 2)

legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth",
                  "Quarter period",
                  "Half period",
                  "Full period"),
       lty = 1, cex = 2.2, bty = "n")

### (2) Trace plot
plot.ts(chain_quarter, col = "orange", lwd = 2,
        ylim = range(c(chain_quarter, chain_half, chain_full)),
        axes = FALSE, ann = FALSE)
lines(chain_half, col = "red",  lwd = 2)
lines(chain_full, col = "blue", lwd = 2)

axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()

mtext(expression(Iterations), side = 1, line = 2.8, cex = cex_lab)
mtext("Trace plot",          side = 2, line = 3.5, cex = cex_lab)

legend("topright",
       col = c("orange", "red", "blue"),
       legend = c("Quarter period",
                  "Half period",
                  "Full period"),
       lty = 1, cex = 2, bty = "n")

### (3) ACF plot
lagmax <- 50
acf_full    <- acf(chain_full,    plot = FALSE, lag.max = lagmax)
acf_half    <- acf(chain_half,    plot = FALSE, lag.max = lagmax)
acf_quarter <- acf(chain_quarter, plot = FALSE, lag.max = lagmax)

lags <- as.numeric(acf_full$lag)

plot(lags, as.numeric(acf_full$acf), type = "l", lwd = 2, col = "blue",
     ylim = range(c(acf_full$acf, acf_half$acf, acf_quarter$acf)),
     axes = FALSE, xlab = "", ylab = "")

lines(lags, as.numeric(acf_half$acf),    col = "red",    lwd = 2)
lines(lags, as.numeric(acf_quarter$acf), col = "orange", lwd = 2)

axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()

mtext("Lag", side = 1, line = 2.8, cex = cex_lab)
mtext("ACF", side = 2, line = 3.5, cex = cex_lab)

legend("topright",
       col = c("orange", "red", "blue"),
       legend = c("Quarter period",
                  "Half period",
                  "Full period"),
       lty = 1, cex = 2, bty = "n")

dev.off()


########## Leapfrog trajectories for different L (fixed epsilon) ##########

pdf("leapfrog_traj_periodicity.pdf", width = 16, height = 5)

par(mfrow = c(1, 3), mar = c(4, 5, 1, 1), oma = c(0,0,0,0))

# Plotting settings
cex_lab  <- 1.5
cex_axis <- 1.2
cex_leg  <- 2
pt_cex   <- 1.5
lwd_line <- 2

# True Hamiltonian circular orbit
r <- 3
a <- 1
t_cont <- seq(0, 2 * pi, length.out = 1000)
x_true <- r * cos(a + t_cont)
p_true <- -r * sin(a + t_cont)

# Fixed epsilon
eps <- 0.1

# Three choices of L: quarter, half, full period
L_vals <- c(
  round(pi / (2 * eps)),   # Quarter period
  round(pi / eps),         # Half period
  round(2 * pi / eps)      # Full period
)

cols <- c("darkgreen", "darkblue", "darkred")

for (k in seq_along(L_vals)) {
  
  L <- L_vals[k]
  
  # Leapfrog trajectory arrays
  e_q <- numeric(L)
  e_p <- numeric(L)
  
  # Initial condition matches analytic start
  e_q[1] <- x_true[1]
  e_p[1] <- p_true[1] - (eps/2) * e_q[1]
  
  # Leapfrog iterations
  if (L >= 2) {
    for (i in 2:L) {
      e_q[i] <- e_q[i-1] + eps * e_p[i-1]
      if (i != L) e_p[i] <- e_p[i-1] - eps * e_q[i]
    }
    e_p[L] <- e_p[L-1] - (eps/2) * e_q[L]
  }
  
  # Panel: true circular orbit + leapfrog trajectory
  plot(x_true, p_true, type = "l", lwd = lwd_line, col = "black",
       xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), asp = 1,
       xlab = expression(Position~(x)),
       ylab = expression(Momentum~(p)),
       cex.lab = cex_lab, cex.axis = cex_axis, las = 1)
  
  lines(e_q, e_p, type = "b", pch = 16, cex = pt_cex,
        col = cols[k], lty = 2, lwd = lwd_line)
  
  # Legend: true Hamiltonian + L
  legend("topright",
         legend = c(
           expression(H_true),
           bquote(L == .(L))
         ),
         col = c("black", cols[k]),
         lty = c(1,2),
         pch = c(NA, 16),
         lwd = c(lwd_line, lwd_line),
         bty = "n", cex = cex_leg)
}

dev.off()