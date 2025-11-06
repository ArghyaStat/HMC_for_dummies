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

# ---- Create one list per column of samples_hmc ----
# Each element will contain a vector of samples for one parameter
hmc_lists <- lapply(1:d, function(j) samples_lastk[, j])
names(hmc_lists) <- var_names
# ---- Convert each list element into an Smcmc object ----
Smcmc_list1 <- Smcmc(samples_lastk[,1])

# ---- Font sizes ----
cex_lab  <- 1.6
cex_axis <- 1.4
lwd_line <- 2
lag.max = 1e3

pdf("HMC_plots.pdf", width = 12, height = 5)
par(mfrow = c(1,3))
SimTools::traceplot(Smcmc_list$chains$Intercept, 
                    ylab = "Intercept", 
                    col = "darkorange",
                    cex.lab = cex_lab, cex.axis = cex_axis)
abline(h = mean(samples_hmc[, 1]), col = "firebrick", lwd = 2, lty = 2)
SimTools::acfplot(Smcmc_list1, 
                  lwd = lwd_line, 
                  chain.col = "red",
                  xlab = "Lag",
                  ylab = " ",
                  lag.max = lag.max)
SimTools::densityplot(Smcmc_list1,
                      Q = c(0.05, 0.95))
dev.off()

Smcmc_list$b.size

plot(density(Smcmc_list1$stacked), xlab = "value", main = " ",
     ylab = "Density")
CIs <- getCI(Smcmc_list1, Q = c(0.05, 0.95))
addCI(Smcmc_list1, CIs, component = 1)


SimTools::ACF(Smcmc_list1,component = NULL, type = c("correlation", "covariance"),
              plot= TRUE, lag.max = NULL, avg.col = "blue", chain.col = "red",
              na.action = na.fail, auto.layout = TRUE, ask = dev.interactive())
