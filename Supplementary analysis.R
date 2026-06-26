##############################################################
####   Supplementary analysis for offset article          ####
####   These are generate non-essential data and figures  ####
####   J.A.Smith   NSW DPI   June 2026                    ####
##############################################################

## It produces 3 supplementary outputs:
##   AIC and deviance explained for M1-M6 (deviance explained not very insightful)
##   Poisson vs negative binomial MAE
##   collinearity-gradient sensitivity
##
## Outputs written to ./supplementary_outputs/

library(mgcv)
library(dplyr)
library(tidyr)

## Load the original and supplementary functions
source("Functions.R")
source("Functions_supplementary.R")

out_dir <- "supplementary_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## 1) AIC and deviance explained

aic_tab <- aic_devexpl_table(distribution = "nbinom")
devexpl_wide <- aic_tab %>%
  select(dataset, model, dev_expl) %>%
  pivot_wider(names_from = model, values_from = dev_expl) %>%
  arrange(dataset)
dAIC_wide <- aic_tab %>%
  select(dataset, model, dAIC) %>%
  pivot_wider(names_from = model, values_from = dAIC) %>%
  arrange(dataset)

write.csv(devexpl_wide, file.path(out_dir, "TableS_devianceExplained_wide.csv"),
          row.names = FALSE)
write.csv(dAIC_wide, file.path(out_dir, "TableS_dAIC_wide.csv"), row.names = FALSE)

print(dAIC_wide)


## 2) Poisson vs negative binomial data generation
## This shows that percentage differences among models increase under Poisson
## data (higher signal-to-noise)
dist_cmp <- distribution_cv_comparison(k = 5, n_repeats = 10, seed = 117,
                                       distributions = c("nbinom", "poisson"))

# Combine into one table
dist_tab <- do.call(rbind, lapply(names(dist_cmp), function(d) {
  m <- as.data.frame(dist_cmp[[d]])
  m$distribution <- d
  m$effort_group <- rownames(m)
  m
}))
rownames(dist_tab) <- NULL
write.csv(dist_tab, file.path(out_dir, "TableS_distribution_MAE_pctDiff.csv"),
          row.names = FALSE)

# Negatuve binomial
print(dist_cmp$nbinom)
# Poisson
print(dist_cmp$poisson)

## 3) Collinearity-gradient sensitivity analysis (Reviewer 2)
## Varies effort-temperature covariance in a single-collinearity proportional
## scenario
coll <- collinearity_gradient(cov_seq = seq(0, 0.9, by = 0.1),
                              n_obs = 1200, k = 5, n_repeats = 10, seed = 117)
write.csv(coll, file.path(out_dir, "TableS_collinearity_gradient.csv"),
          row.names = FALSE)
print(coll)

## 3b) Better to replicate the above (de-noises the slope recovery)
## This repeats each covariance value over many datasets. Can be slow
n_rep_coll <- 100  # datasets per covariance value
coll_rep <- collinearity_gradient_replicated(cov_seq = seq(0, 0.9, by = 0.1),
                                             n_rep = n_rep_coll,
                                             n_obs = 1200, base_seed = 117)

write.csv(coll_rep$summary,
          file.path(out_dir, "TableS_collinearity_gradient_replicated.csv"),
          row.names = FALSE)

# Figures: a) M2 slope mean with 2.5-97.5% replicate band, true = 1;
# b) proportion of replicates recovering proportionality.
s <- coll_rep$summary
png(file.path(out_dir, "FigS_collinearity_gradient_replicated.png"),
    width = 1700, height = 750, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))

# a) slope: mean + 2.5-97.5% band across replicates
plot(s$cov, s$slope_mean, type = "n",
     ylim = range(c(s$slope_lwr, s$slope_upr, 1)), las = 1,
     xlab = "Effort-temperature covariance",
     ylab = "Estimated log(effort) slope (M2)",
     main = sprintf("(a) Slope recovery over %d replicates", n_rep_coll))
polygon(c(s$cov, rev(s$cov)), c(s$slope_lwr, rev(s$slope_upr)),
        col = adjustcolor("blue", 0.15), border = NA)
lines(s$cov, s$slope_mean, type = "b", pch = 16, col = "blue")
abline(h = 1, lty = 2, col = "grey40")  # true proportional value
legend("topright", bty = "n",
       legend = c("mean slope", "2.5-97.5% of replicates", "true value (slope = 1)"),
       col = c("blue", adjustcolor("blue", 0.4), "grey40"),
       pch = c(16, 15, NA), lty = c(1, NA, 2))

# b) proportion of replicates whose 95% CI includes the true slope of 1
plot(s$cov, s$prop_recover, type = "b", pch = 16, col = "darkgreen",
     ylim = c(0, 1), las = 1,
     xlab = "Effort-temperature covariance",
     ylab = "Proportion of replicates recovering slope = 1",
     main = "(b) Recovery rate")
abline(h = 0.95, lty = 3, col = "grey50")  # nominal 95% coverage
legend("bottomleft", bty = "n", lty = 3, col = "grey50",
       legend = "nominal 0.95 coverage")

dev.off()
