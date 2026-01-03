rm(list = ls())

mydir <- this.path::here()
setwd(mydir)
set.seed(1)

# Figure 1 in Section 3
pdf("gaussian_potential.pdf", width = 12, height = 5)

x <- seq(-3, 3, length = 5e2)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

# Left plot: Density
plot(x, dnorm(x), type = 'l',
     ylab = expression(Density~pi(x)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.2,    # enlarge axis labels
     cex.axis = 1.2)   # enlarge tick labels
points(x = -2, y = 0, pch = 16, col = "blue", cex = 1.3)

# Right plot: Potential
plot(x, -dnorm(x, log = TRUE), type = 'l',
     ylab = expression(Potential~U(x)),
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.2,
     cex.axis = 1.2)
points(x = -2, y = -dnorm(-2, log = TRUE), pch = 16, col = "blue", cex = 1.3)

dev.off()

### Hamiltonian dynamics illustration
# Figure 2 in Section 3

pdf("hamiltnoian_dynamics1.pdf", width = 12, height = 5)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

x0 <- 2
p0 <- -2
t <- seq(0, 2*pi, length = 5e2)
x <- x0*cos(t) + p0*sin(t)
p <- -x0*sin(t) + p0*cos(t)


plot(x, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3),
     asp = 1, lty = 2,  ylab = expression(Momentum~(p)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.3,    # enlarge axis labels
     cex.axis = 1.2, lwd = 2)
lines(x[1:101], p[1:101], lwd = 2)
points(x[1], p[1], col = "blue", pch = 16, cex = 1.3)
text(x[1] + .6, p[1] - .05, "t = 0")

points(x[101], p[101], col = "blue", pch = 16, cex = 1.3)
text(x[101] + .6, p[101] - .05, "t = 2")
qt <- seq(-3, 3, length = 1e3)

plot(x, -dnorm(x, log = TRUE), type = 'l',  
     ylab = expression(Potential~U(x)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.3,    # enlarge axis labels
     cex.axis = 1.2, lwd = 2)
lines(x[1:101], -dnorm(x[1:101], log = TRUE), col = "blue", lwd = 2)
points(x = x[1], y = -dnorm(x[1], log = TRUE), pch = 16, col = "blue", cex = 1.3)
text(x[1] + .6, -dnorm(x[1], log = TRUE), "t = 0")
points(x = x[101], y = -dnorm(x[101], log = TRUE), pch = 16, col = "blue", cex = 1.3)

text(x[101] + .6, -dnorm(x[101], log = TRUE), "t = 2")
dev.off()

# Figure 3 in Section 3
pdf("hamiltnoian_dynamics2.pdf", width = 12, height = 5)

# Adjust margins and font sizes
par(mfrow = c(1,2), mar = c(4,5,1,1), oma = c(0,0,0,0))

x_new <- x
p_new <- p
p_new[1] <- 1

p0 <- -1
t <- seq(0, 2*pi, length = 5e2)
x_new <- x0*cos(t) + p0*sin(t)
p_new <- -x0*sin(t) + p0*cos(t)


plot(x, p, type = 'l', ylim = c(-3,3), xlim = c(-3,3), 
     asp = 1, lty = 2, ylab = expression(Momentum~(p)),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2, lwd = 2
     )
lines(x[1:101], p[1:101], lwd = 2)
points(x[1], p[1], col = "blue", pch = 16, cex = 1.3)

lines(x_new, p_new, col = "purple", lty = 2, lwd = 2)
lines(x_new[1:101], p_new[1:101], col = "purple", lty = 1, lwd = 2)
points(x_new[1], p_new[1], col = "purple", pch = 16, cex = 1.3)
text(x_new[1]-.6, p_new[1] - .05, "t = 0")
points(x[101], p[101], col = "blue", pch = 16, cex = 1.3)

points(x_new[101], p_new[101], col = "purple", pch = 16)
text(x_new[101]+.2, p_new[101] + .2, "t = 2")


qt <- seq(-3, 3, length = 1e3)

plot(x, -dnorm(x, log = TRUE), type = 'l', 
     ylab = expression(Potential~(U(x))),   # ← formula now rendered
     xlab = expression(Position~(x)), las = 1,
     cex.lab = 1.5,    # enlarge axis labels
     cex.axis = 1.2, lwd = 2)
lines(x_new[1:101], -dnorm(x_new[1:101], log = TRUE), col = "purple", lwd = 2)
points(x = x_new[1], y = -dnorm(x_new[1], log = TRUE), pch = 16, col = "purple", cex = 1.3)
text(x_new[1]+.4, -dnorm(x_new[1], log = TRUE), "t = 0")
points(x = x_new[101], y = -dnorm(x_new[101], log = TRUE), pch = 16, col = "purple", cex = 1.3)
text(x_new[101]+.3, -dnorm(x_new[101], log = TRUE) + .1, "t = 2")

dev.off()


##### Gaussian HMC
# idealized HMC function
normal_idealHMC <- function(s = 1, n = 1e4)
{
  x <- numeric(length = n)
  x[1] <- 0 # starting at zero position
  # don't need a starting value of p
  for(k in 2:n)
  {
    p <- rnorm(1)   # new momentum
    
    # simulating Hamiltonian forward s time units 
    # no practical need to flip momentum
    x[k] <- x[k-1]*cos(s) + p*sin(s)
    p <- -x[k-1]*sin(s) + p*cos(s)
  }
  return(x) # not returning momentum
}


chain1 <- normal_idealHMC(s = .1)
chain2 <- normal_idealHMC(s = 1)
chain3 <- normal_idealHMC(s = 5)

# Figure 4 in Section 3
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
lines(density(chain1), col = "blue",   lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain3), col = "orange", lwd = 2)
legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth", "s = 5", "s = 1", "s = .1"),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot: suppress axes/labels and add them manually with cex
plot.ts(chain1, col = "blue", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE,  cex.lab = cex_lab)   # suppress default axes and labels
lines(chain2, col = "red",  lwd = 2)
lines(chain3, col = "orange", lwd = 2)
# add axes with explicit cex.axis
axis(1, cex.axis = cex_axis)   # x-axis
axis(2, cex.axis = cex_axis)   # y-axis
box()
# add axis labels with explicit cex
mtext(expression(Iterations), side = 1, line = 2.8, cex = 1.2)
mtext(expression("Trace plot"), side = 2, line = 3.5, cex = 1.2)
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
plot(lags, as.numeric(acf3$acf), type = "l", lwd = 2, col = "orange",
     ylim = c(min(c(acf1$acf, acf2$acf, acf3$acf)),
              max(c(acf1$acf, acf2$acf, acf3$acf))),
     axes = FALSE, xlab = "", ylab = "", cex.lab = 1.2)
lines(lags, as.numeric(acf2$acf), col = "red",    lwd = 2)
lines(lags, as.numeric(acf1$acf), col = "blue", lwd = 2)

# add axes with explicit sizes
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext("Lag", side = 1, line = 2.8, cex = 1.2)
mtext("ACF", side = 2, line = 3.5, cex = 1.2)

legend("topright",
       col = c("orange", "red", "blue"),
       legend = c("s = 5", "s = 1", "s = .1"),
       lty = 1, cex = 2, bty = "n")

dev.off()

# Figure 5 in Section 3
## Leapfrog integrator — one panel per epsilon
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


########## Leapfrog for standard Normal ############


normal_HMC <- function(L = 10, eps = .1, n = 1e4)
{
  x <- numeric(length = n)
  x[1] <- 0 # starting from x = 0
  accept <- 1
  
  # vectorizing this outside to save time in R
  momentums <- rnorm(n) 
  
  for(k in 2:n)
  {
    p <- momentums[k] # new momentum
    x_prop <- x[k-1]
    
    # one half-Euler step
    p_prop <- p - eps/2 * x_prop  
    # let the frog leap!
    for(l in 1:L)
    {
      x_prop <- x_prop + eps * p_prop
      if(l != L) p_prop <- p_prop - eps*x_prop # exception for the last step
    }
    # one last half-Euler step
    p_prop <- p_prop - eps/2 * x_prop
    
    # Accept-reject
    log.ratio <-  (-x_prop^2 - p_prop^2 + x[k-1]^2 + p^2)/2
    if(log(runif(1)) < log.ratio)
    {
      x[k] <- x_prop
      accept <- accept + 1
    } else{
      x[k] <- x[k-1]
    }
  }
  print(paste("Acceptance = ", accept/n))
  return(x)
}

#Keeping L * epsilon = 1
chain1 <- normal_HMC(L = 10, eps = .1)
chain2 <- normal_HMC(L = 1, eps = 1)
chain3 <- normal_HMC(L = 100, eps = .01) # more time consuming

# Figure 6 in Section 4
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
lines(density(chain1), col = "blue",   lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain3), col = "orange", lwd = 2)
legend("topright",
       col = c("black", "red", "blue", "orange"),
       legend = c("Truth", 
                  expression(epsilon == 1),
                  expression(epsilon == 0.1),
                  expression(epsilon == 0.01)),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot
plot.ts(chain1, col = "blue", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE)
lines(chain2, col = "red",  lwd = 2)
lines(chain3, col = "orange", lwd = 2)
axis(1, cex.axis = cex_axis)
axis(2, cex.axis = cex_axis)
box()
mtext(expression(Iterations), side = 1, line = 2.8, cex = cex_lab)
mtext("Trace plot", side = 2, line = 3.5, cex = cex_lab)
legend("topright",
       col = c("red", "blue", "orange"),
       legend = c(expression(epsilon == 1),
                  expression(epsilon == 0.1),
                  expression(epsilon == 0.01)),
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
       col = c("red", "blue", "orange"),
       legend = c(expression(epsilon == 1),
                  expression(epsilon == 0.1),
                  expression(epsilon == 0.01)),
       lty = 1, cex = 2, bty = "n")

dev.off()

# Keeping L * epsilon = 10
chain1 <- normal_HMC(L = 100, eps = .1)
chain2 <- normal_HMC(L = 10,  eps = 1)
chain3 <- normal_HMC(L = 1,   eps = 10)

# Figure 7 in Section 4 
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
lines(density(chain1), col = "blue",   lwd = 2)
lines(density(chain2), col = "red",    lwd = 2)
lines(density(chain3), col = "orange", lwd = 2)
legend("topright",
       col = c("black", "orange", "red", "blue"),
       legend = c("Truth",
                  expression(epsilon == 10),
                  expression(epsilon == 1),
                  expression(epsilon == 0.1)),
       lty = 1, cex = 2.2, bty = "n")

## (2) Trace plot
plot.ts(chain1, col = "blue", lwd = 2, ylim = c(-4, 6),
        axes = FALSE, ann = FALSE)
lines(chain2, col = "red",  lwd = 2)
lines(chain3, col = "orange", lwd = 2)
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

chain_full    <- normal_HMC(L = L_full,    eps = eps)
chain_half    <- normal_HMC(L = L_half,    eps = eps)
chain_quarter <- normal_HMC(L = L_quarter, eps = eps)


# Figure 8 in Section 4
# Leapfrog trajectories for different L (fixed epsilon) #

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