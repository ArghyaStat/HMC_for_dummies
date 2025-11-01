# ============================================================
# Bayesian Logistic Regression using HMC and MH
# (Pima.tr dataset from MASS)
# ============================================================

rm(list = ls())
set.seed(123)

# ---- Load Data ----
data(Pima.tr, package = "MASS")

# ---- Prepare Response and Design Matrix ----
y <- as.numeric(Pima.tr$type == "Yes")
X <- as.matrix(Pima.tr[, !names(Pima.tr) %in% "type"])
X <- cbind(1, X)  # add intercept
colnames(X)[1] <- "(Intercept)"
n <- nrow(X)
d <- ncol(X)  # number of regression coefficients

# ---- Prior ----
sigma2_prior <- 100  # (sigma_prior = 10)^2

# ---- Potential Energy and Gradient ----
U <- function(beta) {
  eta <- X %*% beta
  loglik <- sum(y * eta - log(1 + exp(eta)))
  prior <- sum(beta^2) / (2 * sigma2_prior)
  return(-loglik + prior)
}

grad_U <- function(beta) {
  eta <- X %*% beta
  p_hat <- 1 / (1 + exp(-eta))
  grad_loglik <- t(X) %*% (y - p_hat)
  grad_prior <- beta / sigma2_prior
  return(-grad_loglik + grad_prior)
}

# ============================================================
# ---- Hamiltonian Monte Carlo (HMC) ----
# ============================================================

HMC <- function(U, grad_U, epsilon, L, current_q) {
  q <- current_q
  p <- rnorm(length(q), 0, 1)  # momentum
  current_p <- p
  
  # Half step for momentum
  p <- p - epsilon * grad_U(q) / 2
  
  # Full steps
  for (i in 1:L) {
    q <- q + epsilon * p
    if (i != L) p <- p - epsilon * grad_U(q)
  }
  
  # Final half step
  p <- p - epsilon * grad_U(q) / 2
  p <- -p  # negate for symmetry
  
  # Hamiltonian energies
  current_U <- U(current_q)
  current_K <- sum(current_p^2) / 2
  proposed_U <- U(q)
  proposed_K <- sum(p^2) / 2
  
  # Accept/reject
  accept_prob <- exp(current_U - proposed_U + current_K - proposed_K)
  accept <- runif(1) < accept_prob
  list(q = if (accept) q else current_q, accept = accept)
}

# ============================================================
# ---- Metropolis–Hastings (MH) ----
# ============================================================

MH <- function(U, current_q, proposal_sd) {
  q_prop <- current_q + rnorm(length(current_q), 0, proposal_sd)
  current_U <- U(current_q)
  proposed_U <- U(q_prop)
  accept_prob <- exp(current_U - proposed_U)
  accept <- runif(1) < accept_prob
  list(q = if (accept) q_prop else current_q, accept = accept)
}

# ============================================================
# ---- Run Both Samplers ----
# ============================================================

# --- HMC parameters ---
n_iter_hmc <- 1e5
epsilon <- 1e-3
L <- 20

# --- MH parameters ---
n_iter_mh <- 1e5
proposal_sd <- 5e-3

# --- Initialize ---
beta_hmc <- rep(0, d)
beta_mh  <- rep(0, d)

samples_hmc <- matrix(NA, n_iter_hmc, d)
samples_mh  <- matrix(NA, n_iter_mh, d)
accept_hmc <- numeric(n_iter_hmc)
accept_mh  <- numeric(n_iter_mh)

# ---- Run HMC ----
for (iter in 1:n_iter_hmc) {
  res <- HMC(U, grad_U, epsilon, L, beta_hmc)
  beta_hmc <- res$q
  accept_hmc[iter] <- res$accept
  samples_hmc[iter, ] <- beta_hmc
}

# ---- Run MH ----
for (iter in 1:n_iter_mh) {
  res <- MH(U, beta_mh, proposal_sd)
  beta_mh <- res$q
  accept_mh[iter] <- res$accept
  samples_mh[iter, ] <- beta_mh
}

# ---- Acceptance Rates ----
acc_rate_hmc <- mean(accept_hmc)
acc_rate_mh  <- mean(accept_mh)
cat("=====================================\n")
cat(sprintf("HMC acceptance rate: %.3f\n", acc_rate_hmc))
cat(sprintf("MH  acceptance rate: %.3f\n", acc_rate_mh))
cat("=====================================\n")

# ---- Save Posterior Samples ----
save(samples_hmc, samples_mh, acc_rate_hmc, acc_rate_mh,
     file = "posterior_samples_HMC_MH.RData")

# ============================================================
# ---- Trace Plots ----
# ============================================================

par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (j in 1:d) {
  plot(samples_hmc[, j], type = "l", col = "blue", lwd = 1.3,
       main = paste("HMC Trace:", colnames(X)[j]),
       xlab = "Iteration", ylab = expression(beta))
}

par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (j in 1:d) {
  plot(samples_mh[, j], type = "l", col = "red", lwd = 1.3,
       main = paste("MH Trace:", colnames(X)[j]),
       xlab = "Iteration", ylab = expression(beta))
}

# ============================================================
# ---- Overlapping Posterior Densities ----
# ============================================================

par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (j in 1:d) {
  dens_hmc <- density(samples_hmc[, j])
  dens_mh  <- density(samples_mh[, j])
  xlim <- range(c(dens_hmc$x, dens_mh$x))
  ylim <- range(c(dens_hmc$y, dens_mh$y))
  
  plot(dens_hmc, col = "blue", lwd = 2, xlim = xlim, ylim = ylim,
       main = colnames(X)[j], xlab = "", ylab = "")
  lines(dens_mh, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("HMC", "MH"),
         col = c("blue", "red"), lwd = 2, lty = c(1, 2),
         bty = "n", cex = 0.8)
}

# ============================================================
# ---- Posterior Means ----
# ============================================================

cat("\nPosterior means (HMC):\n")
print(round(apply(samples_hmc, 2, mean), 3))
cat("\nPosterior means (MH):\n")
print(round(apply(samples_mh, 2, mean), 3))

glm_init <- glm.fit(X, y, family = binomial())
beta_init <- glm_init$coefficients
