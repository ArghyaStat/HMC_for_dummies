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
sigma2_beta <- 1e2  # (sigma_prior = 10)^2

# ---- Potential Energy and Gradient ----
U <- function(beta) {
  eta <- X %*% beta
  loglik <- sum(y * eta - log(1 + exp(eta)))
  prior <- sum(beta^2) / (2 * sigma2_beta)
  return(-loglik + prior)
}

grad_U <- function(beta) {
  eta <- X %*% beta
  p_hat <- 1 / (1 + exp(-eta))
  grad_loglik <- t(X) %*% (y - p_hat)
  grad_prior <- beta / sigma2_beta
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
  
  log_accept_ratio <- (current_U - proposed_U) + (current_K - proposed_K)
  
  # ---- Metropolis-Hastings acceptance step (log-scale) ----
  if (log(runif(1)) < log_accept_ratio) {
    list(q = q, accept = TRUE)
  } else {
    list(q = current_q, accept = FALSE)
  }
}




# --- HMC parameters ---
niters <- 2e5
epsilon <- 1e-3
L <- 20


# --- Initialize beta ---
beta_hmc <- rep(0, d)

samples_hmc <- matrix(NA, niters, d)
accept_hmc <- numeric(niters)


# ---- Run HMC ----
for (iter in 1:niters) {
  res <- HMC(U, grad_U, epsilon, L, beta_hmc)
  beta_hmc <- res$q
  accept_hmc[iter] <- res$accept
  samples_hmc[iter, ] <- beta_hmc
}


# ---- Acceptance Rates ----
acc_rate_hmc <- mean(accept_hmc)
cat(sprintf("HMC acceptance rate: %.3f\n", acc_rate_hmc))

# ---- Save Posterior Samples ----
save(samples_hmc, acc_rate_hmc, file = "posterior_samples_HMC.RData")

##### Posterior Diagnostics for HMC: Density, Trace, and ACF #####

##### HMC Density Plots for Posterior Samples #####

# ---- Load posterior samples ----
load("posterior_samples_HMC.RData")   # loads: samples_hmc, acc_rate_hmc

# ---- Variable Names ----
var_names <- c("Intercept", "pregnant", "glucose", "pressure",
               "triceps", "mass", "pedigree", "age")

# ---- PDF output ----
pdf("HMC_density_plots.pdf", width = 10, height = 12)

# ---- Layout and margins ----
par(mfrow = c(4, 2),                      # 8 panels → 4 rows × 2 columns
    mar  = c(3.5, 4.5, 0.5, 0.5),         # small margins (bottom, left, top, right)
    oma  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2.2, 0.8, 0))

# ---- Font sizes ----
cex_lab  <- 2
cex_axis <- 1.8
lwd_line <- 2

# ---- Plot densities ----
for (j in seq_len(ncol(samples_hmc))) {
  d <- density(samples_hmc[, j])
  
  plot(d$x, d$y, type = "l", lwd = lwd_line, col = "blue",
       xlab = "Value",
       ylab = var_names[j],
       cex.lab = cex_lab, cex.axis = cex_axis)
  
  abline(v = mean(samples_hmc[, j]), col = "firebrick", lwd = 2, lty = 2)
  
  legend("topright",
         legend = c("Posterior mean"),
         col = "firebrick", lty = 2, lwd = 2,
         bty = "n", cex = 1.4)
}

dev.off()

# ---- User-defined number of last iterations ----
k <- 5e4   # change as needed

n_total <- nrow(samples_hmc)
if (k > n_total) {
  warning(sprintf("Requested k=%d exceeds total samples (%d); using all samples.", k, n_total))
  k <- n_total
}

# ---- Subset last k samples ----
samples_lastk <- tail(samples_hmc, k)

# ---- PDF output ----
pdf("HMC_trace_plots.pdf", width = 10, height = 12)

# ---- Layout and margins ----
par(mfrow = c(4, 2),
    mar  = c(3.5, 4.5, 1.2, 0.5),
    oma  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2.2, 0.8, 0))

# ---- Font sizes ----
cex_lab  <- 2
cex_axis <- 1.8
lwd_line <- 1.5

# ---- Plot trace for each beta ----
for (j in seq_len(ncol(samples_lastk))) {
  beta_j <- samples_lastk[, j]
  mean_j <- mean(samples_lastk[, j])
  
  plot(seq_len(k), beta_j, type = "l", lwd = lwd_line,
       col = "darkorange", bty = "l",
       xlab = "Iterations", ylab = var_names[j],
       cex.lab = cex_lab, cex.axis = cex_axis)
  
  abline(h = mean(samples_hmc[, j]), col = "firebrick", lwd = 2, lty = 2)

  legend("topright",
         legend = c("Post.mean"),
         col = "firebrick", lty = 2, lwd = 2,
         bty = "n", cex = 1.4)

  
}

dev.off()

# ---- PDF output ----
pdf("HMC_acf_plots.pdf", width = 10, height = 12)

# ---- Layout and margins ----
par(mfrow = c(4, 2),
    mar  = c(3.5, 4.5, 0.5, 0.5),
    oma  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2.2, 0.8, 0))

# ---- Font sizes ----
cex_lab  <- 2
cex_axis <- 1.8
lwd_line <- 2

# ---- Lag length ----
lag_max <- 1e3

# ---- Plot ACF for each beta ----
for (j in seq_len(ncol(samples_hmc))) {
  acf_obj <- acf(samples_hmc[, j], plot = FALSE, lag.max = lag_max)
  lags <- as.numeric(acf_obj$lag)
  acf_vals <- as.numeric(acf_obj$acf)
  
  plot(lags, acf_vals, type = "l", lwd = lwd_line, col = "darkred",
       xlab = "Lag",
       ylab = var_names[j],
       cex.lab = cex_lab, cex.axis = cex_axis,
       ylim = c(0, 1))
 
}

dev.off()

library(xtable)

if (!exists("samples_hmc")) stop("samples_hmc not found in workspace. Load posterior_samples_HMC.RData first.")
# coerce to matrix
samples_hmc <- as.matrix(samples_hmc)
d <- ncol(samples_hmc)


if (is.null(var_names)) {
  var_names <- if (!is.null(colnames(samples_hmc))) colnames(samples_hmc) else paste0("beta", seq_len(d))
} else {
  # ensure same length as d
  if (length(var_names) == d) {
    var_names <- var_names
  } else if (length(var_names) > d) {
    warning("Provided var_names longer than number of parameters in samples_hmc. Truncating to first ", d, " names.")
    var_names <- var_names[1:d]
  } else { # shorter
    warning("Provided var_names shorter than number of parameters in samples_hmc. Using provided names for first ",
            length(var_names), " and generating generic names for remaining.")
    var_names <- c(var_names, paste0("beta", seq_len(d - length(var_names))))
  }
}

# ---- Compute summaries (safe, base R) ----
means  <- colMeans(samples_hmc)
sds    <- apply(samples_hmc, 2, sd)
lower2 <- apply(samples_hmc, 2, quantile, probs = 0.025, names = FALSE)
upper2 <- apply(samples_hmc, 2, quantile, probs = 0.975, names = FALSE)

# ---- Build tibble/data.frame (no size mismatch) ----
posterior_summary <- data.frame(
  Variable = var_names,
  Mean = means,
  SD = sds,
  Lower = lower2,
  Upper = upper2,
  stringsAsFactors = FALSE
)

# ---- Format numeric columns for LaTeX (keep plain numeric columns too) ----
# If you want \hphantom{0} alignment in the LaTeX output, we create a formatted text version.
fmt_num_phantom0 <- function(x, digits = 3) {
  sapply(x, function(v) {
    if (is.na(v)) return("")
    if (v >= 0) sprintf("\\hphantom{0}%.*f", digits, v) else sprintf("%.*f", digits, v)
  }, USE.NAMES = FALSE)
}

# Formatted strings wrapped in $...$
posterior_summary$Mean_fmt <- paste0("$", fmt_num_phantom0(posterior_summary$Mean, digits = 3), "$")
posterior_summary$SD_fmt   <- paste0("$", fmt_num_phantom0(posterior_summary$SD, digits = 3), "$")
posterior_summary$CI_fmt   <- paste0(
  "$[",
  fmt_num_phantom0(posterior_summary$Lower, digits = 3),
  ",\\, ",
  fmt_num_phantom0(posterior_summary$Upper, digits = 3),
  "]$"
)

# ---- Final table to pass to xtable (keep formatted columns as character) ----
tex_table <- data.frame(
  Variable = posterior_summary$Variable,
  `Posterior mean` = posterior_summary$Mean_fmt,
  `Posterior sd`   = posterior_summary$SD_fmt,
  `Posterior $95\\%$ CI` = posterior_summary$CI_fmt,
  stringsAsFactors = FALSE
)

# ---- Create and print xtable (booktabs = FALSE as requested) ----
xt <- xtable(tex_table,
             caption = "Posterior summary (mean, SD, and 95\\% credible intervals) from HMC samples",
             label = "tab:posterior_summary",
             align = c("l", "l", "r", "r", "l"))

print(xt,
      include.rownames = FALSE,
      sanitize.text.function = identity, # keep LaTeX markup
      comment = FALSE,
      booktabs = FALSE)



# library(dplyr)
# library(xtable)
# 
# # ---- Posterior summary table ----
# posterior_summary <- as_tibble(samples_hmc) %>%
#   summarise(
#     Mean  = apply(., 2, mean),
#     SD    = apply(., 2, sd),
#     lower = apply(., 2, quantile, probs = 0.025),
#     upper = apply(., 2, quantile, probs = 0.975)
#   )
# 
# # ---- Combine with variable names and format credible intervals ----
# posterior_summary <- tibble(
#   Variable = colnames(samples_hmc),
#   Mean = as.numeric(posterior_summary$Mean),
#   SD = as.numeric(posterior_summary$SD),
#   `95% CI` = paste0("[", 
#                     formatC(posterior_summary$lower, digits = 3, format = "f"), ", ",
#                     formatC(posterior_summary$upper, digits = 3, format = "f"), "]")
# )
# 
# # ---- Create xtable object ----
# xtable_summary <- xtable(
#   posterior_summary,
#   digits = 3,  # single scalar fixes the mismatch issue
#   caption = "Posterior summary (mean, SD, and 95\\% credible intervals) from HMC samples",
#   label = "tab:posterior_summary"
# )
# 
# print(
#   xtable_summary,
#   include.rownames = FALSE,
#   sanitize.text.function = identity,
#   comment = FALSE,
#   booktabs = FALSE  # <- use normal LaTeX rules, no \toprule/\midrule
# )
# 
# # Robust posterior-summary -> LaTeX pipeline
# # Assumes you have saved samples_hmc and acc_rate_hmc in posterior_samples_HMC.RData
# # and that samples_hmc is an (iter x d) numeric matrix or data.frame.
# 
# 
