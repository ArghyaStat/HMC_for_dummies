rm(list = ls())

mydir <- this.path::here()
setwd(mydir)

library(SimTools)
library(MASS)

load(file = "posterior_samples_HMC.RData")

# ---- Variable Names ----
var_names <- c("Intercept", "pregnant", "glucose", "pressure",
               "triceps", "mass", "pedigree", "age")
d <- length(var_names)
niters <- nrow(samples_hmc)

k <- 5e4
# ---- Subset last k samples ----
samples_lastk <- tail(samples_hmc, k)

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
