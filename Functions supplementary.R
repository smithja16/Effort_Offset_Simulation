##############################################################
####   Supplementary analysis functions                   ####
####   Drop-in alongside Functions.R / Analysis.R         ####
####   J.A.Smith  NSW DPI  June 2026                      ####
##############################################################

## Source this file after Functions.R, bc it reuses some functions
## This script provides 3 supplementary analyses:
##   1) AIC and deviance explained for M1-M6
##   2) Poisson vs negative binomial data generation
##   3) Collinearity-gradient sensitivity analysis


## Scenario definitions and helpers

# The 8 scenarios, identical to Analysis.R
supp_scenarios <- data.frame(
  temp_effort_cov = c(0, 0.5, 0, 0.5, 0.5, 0.5, 0, 0.5),
  site_effort_cov = c(0, 0,   0, 0,   0,   -0.5, 0, -0.5),
  effort_effect   = c("proportional", "proportional", "threshold", "threshold",
                      "proportional", "proportional", "constant",  "threshold"),
  temp_effect     = c("linear", "linear", "linear", "linear",
                      "nonlinear", "linear", "linear", "nonlinear")
)

# Generate one data scenario (uses the same seed = 117 as Analysis.R)
make_scenario_data <- function(dd, n_obs = 1200, distribution = "nbinom", seed = 117) {
  generate_counts(n_obs = n_obs,
                  temp_effort_cor = supp_scenarios$temp_effort_cov[dd],
                  site_effort_cor = supp_scenarios$site_effort_cov[dd],
                  effort_type     = supp_scenarios$effort_effect[dd],
                  temp_effect     = supp_scenarios$temp_effect[dd],
                  temp_optimum    = 20,
                  temp_breadth    = 5,
                  beta_temp       = 0.1,
                  distribution    = distribution,
                  seed            = seed)
}

# Fit the six models to a data frame; identical specification to Analysis.R, 
# except that we use "ML" to better estimate AIC
fit_models_data <- function(datax) {
  datax$logEffort <- log(datax$Effort)
  list(
    M1 = mgcv::gam(Count ~ Temperature + Site + offset(logEffort),
                   family = mgcv::nb(), data = datax, method="ML"),
    M2 = mgcv::gam(Count ~ Temperature + Site + logEffort,
                   family = mgcv::nb(), data = datax, method="ML"),
    M3 = mgcv::gam(Count ~ Temperature + Site + Effort,
                   family = mgcv::nb(), data = datax, method="ML"),
    M4 = mgcv::gam(Count ~ Temperature + Site + s(logEffort),
                   family = mgcv::nb(), data = datax, method="ML"),
    M5 = mgcv::gam(Count ~ Temperature + Site + s(Effort),
                   family = mgcv::nb(), data = datax, method="ML"),
    M6 = mgcv::gam(Count ~ s(Temperature) + Site + s(logEffort),
                   family = mgcv::nb(), data = datax, method="ML")
  )
}

# Helper: summarise a within-dataset percent-difference-from-minimum table,
# averaged over the proportional (D1,D2,D5,D6) and threshold/constant
# (D3,D4,D7,D8) groups, excluding M6 (matches out_sample_mae in Functions.R).
pct_diff_summary <- function(mae_wide) {
  # mae_wide: data frame with column 'dataset' and one column per model (M1..M5)
  mod_cols <- setdiff(names(mae_wide), "dataset")
  mae_wide[, mod_cols] <- t(apply(mae_wide[, mod_cols], 1, function(x) {
    100 * (x - min(x)) / min(x)
  }))
  prop   <- colMeans(mae_wide[mae_wide$dataset %in% c(1, 2, 5, 6), mod_cols, drop = FALSE])
  thresh <- colMeans(mae_wide[mae_wide$dataset %in% c(3, 4, 7, 8), mod_cols, drop = FALSE])
  out <- round(rbind(Proportional = prop, `Threshold/constant` = thresh), 2)
  return(out)
}


## 1) AIC and deviance explained
## M1-M6 share the same response (Count), family (NB) and data, so AIC and
## deviance explained are directly comparable within a data scenario.
## Returns a long table and writes a wide delta-AIC table

aic_devexpl_table <- function(distribution = "nbinom") {
  res <- expand.grid(dataset = 1:8, model = paste0("M", 1:6),
                     stringsAsFactors = FALSE)
  res$AIC <- NA_real_
  res$dev_expl <- NA_real_
  for (dd in 1:8) {
    datax <- make_scenario_data(dd, distribution = distribution)
    mods  <- fit_models_data(datax)
    for (mm in 1:6) {
      Mx <- mods[[paste0("M", mm)]]
      r  <- which(res$dataset == dd & res$model == paste0("M", mm))
      res$AIC[r]      <- AIC(Mx)
      res$dev_expl[r] <- summary(Mx)$dev.expl
    }
  }
  # delta AIC within each data scenario (lower AIC = better)
  res <- res %>%
    dplyr::group_by(dataset) %>%
    dplyr::mutate(dAIC = round(AIC - min(AIC), 1)) %>%
    dplyr::ungroup()
  res$dev_expl <- round(100 * res$dev_expl, 1)  # as a percentage
  return(as.data.frame(res))
}


## 2) Poisson vs negative binomial data generation
## Re-runs the cross-validation MAE comparison (as in Table 3) under each
## data-generating distribution. The fitted models are unchanged (NB()),
## so only the signal-to-noise of the data changes.

distribution_cv_comparison <- function(k = 5, n_repeats = 10, seed = 117,
                                        distributions = c("nbinom", "poisson")) {
  results <- list()
  for (dist in distributions) {
    overall <- list()
    for (dd in 1:8) {
      datax <- make_scenario_data(dd, distribution = dist)
      datax$logEffort <- log(datax$Effort)
      mods <- fit_models_data(datax)
      cv <- cv_repeated_kfold(mods, datax, k = k, n_repeats = n_repeats, seed = seed)
      ov <- cv$overall
      ov$dataset <- dd
      overall[[dd]] <- ov
    }
    overall <- do.call(rbind, overall)

    # wide table of overall MAE (exclude M6, as in out_sample_mae)
    mae_wide <- overall[overall$model != "M6", c("model", "dataset", "mae_mean_overall")]
    mae_wide <- mae_wide %>%
      tidyr::pivot_wider(names_from = model, values_from = mae_mean_overall) %>%
      as.data.frame()
    results[[dist]] <- pct_diff_summary(mae_wide)
  }
  return(results)
}


## 3) Collinearity-gradient sensitivity analysis

## Varies the effort-temperature covariance in a single-collinearity,
## proportional, linear-temperature scenario (like D2). Only one covariance
## is varied, so the correlation matrix stays positive-definite up to ~0.95.
## For each value it records the CV MAE of the offset (M1), log-covariate (M2)
## and smoother (M4), and the M2 log(effort) slope estimate (true value = 1).

collinearity_gradient <- function(cov_seq = seq(0, 0.9, by = 0.1),
                                   n_obs = 1200, k = 5, n_repeats = 10, seed = 117) {
  out <- data.frame()
  for (cv in cov_seq) {
    datax <- generate_counts(n_obs = n_obs,
                             temp_effort_cor = cv,
                             site_effort_cor = 0,
                             effort_type = "proportional",
                             temp_effect = "linear",
                             beta_temp = 0.1,
                             distribution = "nbinom",
                             seed = seed)
    datax$logEffort <- log(datax$Effort)
    mods <- fit_models_data(datax)

    # cross-validated MAE for the three focal models
    cvres <- cv_repeated_kfold(mods[c("M1", "M2", "M4")], datax,
                               k = k, n_repeats = n_repeats, seed = seed)
    mae <- setNames(cvres$overall$mae_mean_overall, cvres$overall$model)

    # M2 log(effort) slope and 95% CI (true value = 1)
    # (extracted as in save_plot_results: coefficients + summary()$se)
    est <- mods$M2$coefficients["logEffort"]
    se  <- summary(mods$M2)$se["logEffort"]
    lwr <- est - 1.96 * se
    upr <- est + 1.96 * se

    out <- rbind(out, data.frame(
      cov        = cv,
      MAE_M1     = unname(mae["M1"]),
      MAE_M2     = unname(mae["M2"]),
      MAE_M4     = unname(mae["M4"]),
      pct_M2_vs_M1 = round(100 * (mae["M2"] - mae["M1"]) / mae["M1"], 2),
      M2_logEffort = round(unname(est), 3),
      M2_CI_lwr  = round(unname(lwr), 3),
      M2_CI_upr  = round(unname(upr), 3),
      recovers_proportional = as.integer(lwr <= 1 & upr >= 1)
    ))
  }
  rownames(out) <- NULL
  return(out)
}


## 3b) Replicate-averaged collinearity gradient

## The single-realization collinearity_gradient() above is noisy, so this
## version repeats each covariance value over many independent datasets
## and summarises the distribution of the M2 log(effort) slope and of the
## offset-vs-covariate predictive difference.
## To stay fast over many replicates it uses in-sample MAE rather than
## cross-validation; the relative M2-vs-M1 comparison is what matters and it
## was already ~0 under CV.
## Per covariance value it reports: mean / SD / 2.5-97.5% range of the slope
## (true value = 1), the proportion of replicates whose 95% CI includes 1
## (this should sit near 0.95 when estimation is unbiased and fall as
## collinearity biases the estimate), and the mean +/- SD percentage
## MAE difference of the covariate (M2) relative to the offset (M1).

collinearity_gradient_replicated <- function(cov_seq = seq(0, 0.9, by = 0.1),
                                             n_rep = 100, n_obs = 1200,
                                             base_seed = 117) {
  per_rep <- data.frame()
  ci <- 0
  for (cv in cov_seq) {
    ci <- ci + 1
    for (rr in seq_len(n_rep)) {
      seed_use <- base_seed + (ci - 1) * n_rep + rr  # unique, reproducible per (cov, rep)
      datax <- generate_counts(n_obs = n_obs,
                               temp_effort_cor = cv,
                               site_effort_cor = 0,
                               effort_type = "proportional",
                               temp_effect = "linear",
                               beta_temp = 0.1,
                               distribution = "nbinom",
                               seed = seed_use)
      datax$logEffort <- log(datax$Effort)

      M1 <- mgcv::gam(Count ~ Temperature + Site + offset(logEffort),
                      family = mgcv::nb(), data = datax)
      M2 <- mgcv::gam(Count ~ Temperature + Site + logEffort,
                      family = mgcv::nb(), data = datax)

      est <- M2$coefficients["logEffort"]
      se  <- summary(M2)$se["logEffort"]
      lwr <- est - 1.96 * se
      upr <- est + 1.96 * se

      mae1 <- MAE(M1$fitted.values, datax$Count)
      mae2 <- MAE(M2$fitted.values, datax$Count)

      per_rep <- rbind(per_rep, data.frame(
        cov          = cv,
        rep          = rr,
        M2_logEffort = unname(est),
        recovers     = as.integer(lwr <= 1 & upr >= 1),
        pct_M2_vs_M1 = 100 * (mae2 - mae1) / mae1
      ))
    }
  }

  summary_tab <- per_rep %>%
    dplyr::group_by(cov) %>%
    dplyr::summarise(
      n_rep        = dplyr::n(),
      slope_mean   = mean(M2_logEffort),
      slope_sd     = sd(M2_logEffort),
      slope_lwr    = quantile(M2_logEffort, 0.025),
      slope_upr    = quantile(M2_logEffort, 0.975),
      prop_recover = mean(recovers),
      pct_mean     = mean(pct_M2_vs_M1),
      pct_sd       = sd(pct_M2_vs_M1),
      .groups = "drop"
    ) %>%
    as.data.frame()

  # round for presentation (keep per_rep raw)
  num <- setdiff(names(summary_tab), c("cov", "n_rep"))
  summary_tab[num] <- round(summary_tab[num], 3)

  return(list(per_rep = per_rep, summary = summary_tab))
}
