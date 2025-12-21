rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

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
      } 
    }
    
    p_prop <- p_prop - (epsilon / 2) * g
    
    
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
L <- 20
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

# ---- Save Posterior Samples ----
save(beta_warmup, accept_warmup, file = "posterior_beta_warmup.RData")

load("posterior_beta_warmup.RData") 

# --- HMC parameters ---
niters <- 1e5
#epsilon <- 1.1 # tuning for full pre-conditioning
epsilon <- 0.11 # tuning for diagonal pre-conditioning
L <- 20
M <- diag(1/diag(cov(beta_warmup))) # diagonal preconditiong
# solve(cov(beta_warmup))  ## full preconditioning


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

# beta_samples <- beta_warmup
# niters <- warmup

var_names <- c("Intercept", "pregnant", "glucose", "pressure",
               "triceps", "mass", "pedigree", "age")

library(SimTools)
library(MASS)

k <- 1e4
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
  
  plot((niters - k + 1):niters, beta_j, type = "l", lwd = lwd_line,
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

pdf("HMC_acf_combined.pdf", width = 12, height = 5)

par(mfrow = c(1,2), mar = c(5,5,1,1), oma = c(0,0,0,0))

lag_max <- 100
lags <- 0:lag_max
p <- ncol(beta_samples)


acf_naive <- acf(beta_warmup[,1], plot = FALSE, lag.max = lag_max)$acf
ylim_all <- range(acf_naive)  

for (j in seq_len(p)) {
  ylim_all <- range(c(
    ylim_all,
    acf(beta_warmup[,j],  plot = FALSE, lag.max = lag_max)$acf,
    acf(beta_samples[,j], plot = FALSE, lag.max = lag_max)$acf
  ))
}


plot(lags, acf_naive, type = "l", lwd = 2, col = "darkblue", lty = 2,
     xlab = "Lag", ylab = "ACF", las = 1,
     cex.lab = 1.5, cex.axis = 1.2,
     ylim = ylim_all)

for (j in seq_len(p)) {
  lines(lags,
        acf(beta_warmup[,j], plot = FALSE, lag.max = lag_max)$acf,
        col = "darkblue", lwd = 2, lty = 1)
}


abline(h = 0, lty = 2)


acf_pre <- acf(beta_samples[,1], plot = FALSE, lag.max = lag_max)$acf

plot(lags, acf_pre, type = "l", lwd = 2, col = "darkred", lty = 1,
     xlab = "Lag", ylab = "ACF", las = 1,
     cex.lab = 1.5, cex.axis = 1.2,
     ylim = ylim_all)

for (j in seq_len(p)) {
  lines(lags,
        acf(beta_samples[,j], plot = FALSE, lag.max = lag_max)$acf,
        col = "darkred", lwd = 2, lty = 1)
}


abline(h = 0, lty = 2)

dev.off()

########## Posterior summary table from pre-conditioned sampler ##########

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
      sanitize.text.function = identity, 
      comment = FALSE,
      booktabs = FALSE)



# legend("topright",
#        legend = c("naive sampler"),
#        col = "darkblue", lty = 1, lwd = 2, bty = "n", cex = 1.2)
# 

# legend("topright",
#        legend = c("preconditioned sampler"),
#        col = "darkred", lty = 1, lwd = 2, bty = "n", cex = 1.2)
# 
# ---- Single combined ACF plot for all components (8 + 8 = 16) ----
# pdf("HMC_acf_combined.pdf", width = 6, height = 5)
# 
# 
# par(mfrow = c(1, 1),
#     mar  = c(4.0, 5.0, 1.0, 0.5),
#     oma  = c(0, 0, 0, 0),
#     mgp  = c(2.2, 0.8, 0))
# 
# # font / line settings
# cex_lab  <- 1.4
# cex_axis <- 1.1
# lwd_line <- 2
# 
# # lag
# lag_max <- 100
# lags <- 0:lag_max
# 
# # number of components
# p <- ncol(beta_samples)
# 
# # allocate storage for ACF values (each column: one component, rows correspond to lags 0:lag_max)
# acf_pre_mat  <- matrix(NA_real_, nrow = lag_max + 1, ncol = p)  # preconditioned (post-warmup) assumed in beta_samples
# acf_naive_mat <- matrix(NA_real_, nrow = lag_max + 1, ncol = p) # naive / warmup assumed in beta_warmup
# 
# # compute ACFs (same lag grid)
# for (j in seq_len(p)) {
#   acf_pre  <- acf(beta_samples[, j],  plot = FALSE, lag.max = lag_max)$acf
#   acf_naive <- acf(beta_warmup[, j], plot = FALSE, lag.max = lag_max)$acf
#   
#   # store: acf(...) returns an array, coerce to numeric vector
#   acf_pre_mat[, j]   <- as.numeric(acf_pre)
#   acf_naive_mat[, j] <- as.numeric(acf_naive)
# }
# 
# # global y-limits across all 16 series
# ylim_all <- range(c(acf_pre_mat, acf_naive_mat), na.rm = TRUE)
# 
# # base plot (use first preconditioned component to initialise the axes)
# plot(lags, acf_pre_mat[, 1], type = "l", lwd = lwd_line, col = "darkred",
#      xlab = "Lag", ylab = "Autocorrelation (ACF)",
#      cex.lab = cex_lab, cex.axis = cex_axis,
#      ylim = ylim_all, xlim = range(lags))
# 
# # add the remaining preconditioned (post-warmup) ACFs (solid darkred)
# if (p > 1) {
#   for (j in 2:p) {
#     lines(lags, acf_pre_mat[, j], lwd = lwd_line, col = "darkred", lty = 1)
#   }
# }
# 
# # add all naive / warmup ACFs (dashed darkblue, lty = 2)
# for (j in seq_len(p)) {
#   lines(lags, acf_naive_mat[, j], lwd = lwd_line, col = "darkblue", lty = 2)
# }
# 
# # reference zero line
# abline(h = 0, lty = 3)
# 
# # legend
# legend("topright",
#        legend = c("preconditioned sampler", "naive sampler"),
#        col    = c("darkred", "darkblue"),
#        lwd    = c(lwd_line, lwd_line),
#        lty    = c(1, 2),
#        bty    = "n",
#        cex    = 1.4)
# 
# dev.off()

# 
# # ---- PDF output ----
# pdf("HMC_acf_combined.pdf", width = 10, height = 12)
# 
# # ---- Layout and margins ----
# par(mfrow = c(4, 2),
#     mar  = c(3.5, 4.5, 0.5, 0.5),
#     oma  = c(0.5, 0.5, 0.5, 0.5),
#     mgp  = c(2.2, 0.8, 0))
# 
# # ---- Font sizes ----
# cex_lab  <- 2
# cex_axis <- 1.8
# lwd_line <- 2
# 
# # ---- Lag length ----
# lag_max <- 1e2
# 
# # ---- Plot ACF for each beta with warmup overlay ----
# for (j in seq_len(ncol(beta_samples))) {
#   
#   acf_main  <- acf(beta_samples[, j], plot = FALSE, lag.max = lag_max)
#   acf_warm  <- acf(beta_warmup[, j], plot = FALSE, lag.max = lag_max)
#   
#   lags      <- as.numeric(acf_main$lag)  # identical for both
#   acf_main_vals <- as.numeric(acf_main$acf)
#   acf_warm_vals <- as.numeric(acf_warm$acf)
#   
#   # combined y-limits
#   ylim_combined <- range(c(acf_main_vals, acf_warm_vals))
#   
#   # main chain ACF
#   plot(lags, acf_main_vals, type = "l", lwd = lwd_line, col = "darkred",
#        xlab = "Lag",
#        ylab = var_names[j],
#        cex.lab = cex_lab, cex.axis = cex_axis,
#        ylim = ylim_combined)
#   
#   # warmup ACF overlay (solid dark blue)
#   lines(lags, acf_warm_vals, lwd = lwd_line, col = "darkblue")
#   
#   abline(h = 0, lty = 3)
#   
#   legend("topright",
#          legend = c("pre-conditioned sampler", "naive sampler"),
#          col = c("darkred", "darkblue"),
#          lwd = lwd_line,
#          lty = 1,
#          bty = "n",
#          cex = 1.4)
# }
# 
# dev.off()

# 
# # ---- PDF output ----
# pdf("HMC_acf_plots.pdf", width = 10, height = 12)
# 
# # ---- Layout and margins ----
# par(mfrow = c(4, 2),
#     mar  = c(3.5, 4.5, 0.5, 0.5),
#     oma  = c(0.5, 0.5, 0.5, 0.5),
#     mgp  = c(2.2, 0.8, 0))
# 
# # ---- Font sizes ----
# cex_lab  <- 2
# cex_axis <- 1.8
# lwd_line <- 2
# 
# # ---- Lag length ----
# lag_max <- 1e2
# 
# # ---- Plot ACF for each beta ----
# for (j in seq_len(ncol(beta_samples))) {
#   acf_obj <- acf(beta_samples[, j], plot = FALSE, lag.max = lag_max)
#   lags <- as.numeric(acf_obj$lag)
#   acf_vals <- as.numeric(acf_obj$acf)
#   
#   plot(lags, acf_vals, type = "l", lwd = lwd_line, col = "darkred",
#        xlab = "Lag",
#        ylab = var_names[j],
#        cex.lab = cex_lab, cex.axis = cex_axis,
#        ylim = c(0, 1))
#   
# }
# 
# dev.off()