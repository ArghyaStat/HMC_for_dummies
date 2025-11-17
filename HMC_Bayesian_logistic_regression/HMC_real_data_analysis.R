rm(list = ls())

# ---- Hamiltonian Monte Carlo (HMC) ----

leapfrog_hmc <- function(epsilon, L, niters, U, grad, x_init, M) {
  
  d <- length(x_init)
  M_inv <- chol2inv(chol(M))
  
  # ---- store samples ----
  x <- matrix(NA, nrow = niters, ncol = d)
  x[1, ] <- x_init
  
  # ---- pre-draw all momentum values ----
  
  chol_M <- chol(M)
  
  p <- matrix(rnorm(niters * d), nrow = niters, ncol = d) %*% chol_M
  
  acc <- 0
  
 
  for (k in 2:niters) {
    
    x_prop <- x[k-1, ]        
    p_prop <- p[k, ]        
    
    g <- grad(x_prop)
    p_prop <- p_prop - (epsilon / 2) * g
    
   
    for (l in 1:L) {
      
      # position update
      x_prop <- x_prop + epsilon * (M_inv %*% p_prop)
      
      # momentum update
      g <- grad(x_prop)
      if (l != L) {
        p_prop <- p_prop - epsilon * g
      } else {
        p_prop <- p_prop - (epsilon / 2) * g
      }
    }
    
   
    x_curr <- x[k-1, ]
    p_curr <- p[k, ]
    
    H_curr <- U(x_curr) + 0.5 * t(p_curr) %*% (M_inv %*% p_curr)
    H_prop <- U(x_prop) + 0.5 * t(p_prop) %*% (M_inv %*% p_prop)
    
    log_alpha <- -(H_prop - H_curr)
    
    # ---- Metropolis accept/reject ----
    if (log(runif(1)) < log_alpha) {
      x[k, ] <- x_prop
      acc <- acc + 1
    } else {
      x[k, ] <- x_curr
    }
  }
  
  list(samples = x, acceptance = acc/(niters))
}

set.seed(123)

# ---- Load Data ----
data(Pima.tr, package = "MASS")

# ---- Prepare Response and Design Matrix ----
y <- as.numeric(Pima.tr$type == "Yes")
X <- as.matrix(Pima.tr[, !names(Pima.tr) %in% "type"])
X <- cbind(1, X)  # add intercept
# colnames(X)[1] <- "Intercept"
n <- nrow(X)
d <- ncol(X)  # number of regression coefficients

# ---- Prior ----
sigma2_beta <- 1e2 

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



# --- HMC parameters ---
warmup <- 1e5
epsilon <- 2.1e-3 
L <- 30
M <- diag(d)


# --- Initialize beta using MLE ---
beta_init <- glm(y ~ X - 1, family = "binomial")$coefficients

hmc_warmup <- leapfrog_hmc(epsilon = epsilon, 
                            L = L, 
                            niters = warmup, 
                            U = U, 
                            grad = grad_U, 
                            x_init = beta_init, 
                            M = M)
accept_warmup <- hmc_warmup$acceptance
cat(sprintf("HMC acceptance rate: %.4f\n", accept_warmup))
beta_warmup <- hmc_warmup$samples



# --- HMC parameters ---
niters <- 1e5
epsilon <- 2.4e-6
L <- 30
M <- cov(beta_warmup)
  # diag(diag(cov(beta_warmup)))


# --- Initialize beta ---
beta_init <- beta_warmup[warmup, ]

hmc_output <- leapfrog_hmc(epsilon = epsilon, 
                           L = L, 
                           niters = niters, 
                           U = U, 
                           grad = grad_U, 
                           x_init = beta_init, 
                           M = M)
accept_main <- hmc_output$acceptance
cat(sprintf("HMC acceptance rate: %.4f\n", accept_main))
beta_samples <- hmc_output$samples


# ---- Save Posterior Samples ----
save(beta_samples, accept_main, file = "posterior_beta_samples.RData")


##### HMC Density Plots for Posterior Samples #####


load("posterior_beta_samples.RData")   

var_names <- c("Intercept", "pregnant", "glucose", "pressure",
               "triceps", "mass", "pedigree", "age")

library(SimTools)
library(MASS)

k <- 1e5
# ---- Subset last k samples ----
samples_lastk <- tail(beta_samples, k)
dim(samples_lastk)

# ---- Create one list per column ----
hmc_lists <- lapply(1:d, function(j) samples_lastk[, j])
names(hmc_lists) <- var_names

# ---- Convert each column to an Smcmc object ----
Smcmc_list <- lapply(hmc_lists, Smcmc)

cex_lab  <- 2
cex_axis <- 1.8

# ---- PDF output ----
pdf("SimTools_density_plots.pdf", width = 10, height = 12)

# ---- Layout and margins ----
par(mfrow = c(4, 2),                      # 8 panels → 4 rows × 2 columns
    mar  = c(3.5, 4.5, 0.5, 0.5),         # margins: bottom, left, top, right
    oma  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2.2, 0.8, 0))                # axis label spacing

# ---- Generate density plots ----
for (j in 1:d) {
  smcmc_obj <- Smcmc_list[[j]]
  dens <- density(smcmc_obj$stacked)
  
  plot(dens,
       main = "",
       xlab = var_names[j],
       ylab = "Density",
       cex.lab = cex_lab, cex.axis = cex_axis)
  
  # Add 90% credible interval
  CIs <- getCI(smcmc_obj, Q = c(0.05, 0.95))
  addCI(smcmc_obj, CIs, component = 1)
}

dev.off()

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
for (j in 1:d) {
  beta_j <- samples_lastk[, j]
  mean_j <- mean(samples_lastk[, j])
  
  plot(1:k, beta_j, type = "l", lwd = lwd_line,
       col = "darkorange", bty = "l",
       xlab = "Iterations", ylab = var_names[j],
       cex.lab = cex_lab, cex.axis = cex_axis)
  
  abline(h = mean(beta_samples[, j]), col = "firebrick", lwd = 2, lty = 2)

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
for (j in seq_len(ncol(beta_samples))) {
  acf_obj <- acf(beta_samples[, j], plot = FALSE, lag.max = lag_max)
  lags <- as.numeric(acf_obj$lag)
  acf_vals <- as.numeric(acf_obj$acf)
  
  plot(lags, acf_vals, type = "l", lwd = lwd_line, col = "darkred",
       xlab = "Lag",
       ylab = var_names[j],
       cex.lab = cex_lab, cex.axis = cex_axis,
       ylim = c(0, 1))
 
}

dev.off()

########## Posterior summary table ##########

library(xtable)

if (!exists("beta_samples")) stop("beta_samples not found in workspace. Load posterior_beta_samples.RData first.")
# coerce to matrix
beta_samples <- as.matrix(beta_samples)
d <- ncol(beta_samples)


if (is.null(var_names)) {
  var_names <- if (!is.null(colnames(beta_samples))) colnames(beta_samples) else paste0("beta", seq_len(d))
} else {
  # ensure same length as d
  if (length(var_names) == d) {
    var_names <- var_names
  } else if (length(var_names) > d) {
    warning("Provided var_names longer than number of parameters in beta_samples. Truncating to first ", d, " names.")
    var_names <- var_names[1:d]
  } else { # shorter
    warning("Provided var_names shorter than number of parameters in beta_samples. Using provided names for first ",
            length(var_names), " and generating generic names for remaining.")
    var_names <- c(var_names, paste0("beta", seq_len(d - length(var_names))))
  }
}

# ---- Compute summaries (safe, base R) ----
means  <- colMeans(beta_samples)
sds    <- apply(beta_samples, 2, sd)
lower2 <- apply(beta_samples, 2, quantile, probs = 0.025, names = FALSE)
upper2 <- apply(beta_samples, 2, quantile, probs = 0.975, names = FALSE)

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
# posterior_summary <- as_tibble(beta_samples) %>%
#   summarise(
#     Mean  = apply(., 2, mean),
#     SD    = apply(., 2, sd),
#     lower = apply(., 2, quantile, probs = 0.025),
#     upper = apply(., 2, quantile, probs = 0.975)
#   )
# 
# # ---- Combine with variable names and format credible intervals ----
# posterior_summary <- tibble(
#   Variable = colnames(beta_samples),
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
# # Assumes you have saved beta_samples and acc_rate_hmc in posterior_beta_samples.RData
# # and that beta_samples is an (iter x d) numeric matrix or data.frame.
# 
# 

# # ---- PDF output ----
# pdf("HMC_density_plots.pdf", width = 10, height = 12)
# 
# # ---- Layout and margins ----
# par(mfrow = c(4, 2),                      # 8 panels → 4 rows × 2 columns
#     mar  = c(3.5, 4.5, 0.5, 0.5),         # small margins (bottom, left, top, right)
#     oma  = c(0.5, 0.5, 0.5, 0.5),
#     mgp  = c(2.2, 0.8, 0))
# 
# # ---- Font sizes ----
# cex_lab  <- 2
# cex_axis <- 1.8
# lwd_line <- 2
# 
# # ---- Plot densities ----
# for (j in 1 : d) {
#   
#   den <- density(beta_samples[, j])
#   
#   density(beta_samples[, j], type = "l", lwd = lwd_line, col = "blue",
#        xlab = "Value",
#        ylab = var_names[j],
#        cex.lab = cex_lab, cex.axis = cex_axis)
#   
#   abline(v = mean(beta_samples[, j]), col = "firebrick", lwd = 2, lty = 2)
#   
#   legend("topright",
#          legend = c("Posterior mean"),
#          col = "firebrick", lty = 2, lwd = 2,
#          bty = "n", cex = 1.4)
# }
# 
# dev.off()

