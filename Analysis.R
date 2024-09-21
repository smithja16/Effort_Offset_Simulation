##############################################################
####   Code to generate and evaluate catch data...        ####
####   ... given Site, Temperature, and Effort effects... ####
####   ... to test how to best model Effort in GLMs.      ####
####   J.A.Smith NSW DPI 30/8/24                          ####
##############################################################


## Load libraries
library(mgcv)
library(dplyr)

## Load functions
source("Functions.R")


## 1) Generate the data scenarios
## Here we generate 8 based on useful changes in collinearity, effort effects, and temp effects

# Number of counts to generate
n_obs <- 1200

# The 8 scenarios as a data frame
scenarios <- data.frame(temp_effort_cov = c(0, 0.5, 0, 0.5, 0.5, 0.5, 0, 0.5),
                        site_effort_cov = c(0, 0, 0, 0, 0, -0.5, 0, -0.5),
                        effort_effect = c("proportional", "proportional", "threshold", "threshold",
                                          "proportional","proportional", "constant", "threshold"),
                        temp_effect = c("linear", "linear", "linear", "linear",
                                          "nonlinear", "linear", "linear", "nonlinear"))

# Create and save the data sets
for (dd in 1:8) {
  counts_dd <- generate_counts(n_obs = n_obs,
                               temp_effort_cor = scenarios$temp_effort_cov[dd],
                               site_effort_cor = scenarios$site_effort_cov[dd],
                               effort_type = scenarios$effort_effect[dd],
                               temp_effect = scenarios$temp_effect[dd],
                               temp_optimum = 20,
                               temp_breadth = 5,
                               beta_temp = 0.1,
                               seed=117)
  assign(paste0("data",dd), counts_dd)
}

# Plot histograms for data scenarios
# note that collinearity can change the magnitude of counts
par(mfrow=c(3,3))
for (dd in 1:8) {
  datax <- get(paste0("data",dd))
  hist(datax$Count, breaks=30,
       main=paste0("Data=",dd))
}


## 2) For each data scenario, fit these six models

fit_models <- function(data_num) {  #e.g. number suffix
  datax <- get(paste0("data",data_num))
  datax$logEffort <- log(datax$Effort)  #pre-transforming makes it easier to get Effort marginal effects
  
  M1 <- gam(Count ~ Temperature + Site + offset(logEffort), family = nb(), data=datax)
  M2 <- gam(Count ~ Temperature + Site + logEffort, family = nb(), data=datax)
  M3 <- gam(Count ~ Temperature + Site + Effort, family = nb(), data=datax)
  M4 <- gam(Count ~ Temperature + Site + s(logEffort), family = nb(), data=datax)
  M5 <- gam(Count ~ Temperature + Site + s(Effort), family = nb(), data=datax)
  M6 <- gam(Count ~ s(Temperature) + Site + s(logEffort), family = nb(), data=datax)

  return(list(M1=M1, M2=M2, M3=M3, M4=M4, M5=M5, M6=M6))
}

for (dd in 1:8) {
  assign(paste0("results_",dd), fit_models(data_num = dd))
}


## 3) Extract results

results <- save_plot_results(ylims=c(200, 200, 150, 150, 60, 200, 15, 50))


## 4) Interpret results

# Evaluate residuals
for (dd in c(1:8)) {
  resultsx <- get(paste0("results_",dd))
  for (mm in 1:6) {
    Mx <- resultsx[[paste0("M",mm)]]
    plot(DHARMa::simulateResiduals(Mx),
         main=paste0("Data=",dd," Model=",mm))
  }
}

# MAE of all data
save_mae <- results$save_mae
wide_mae <- save_mae[save_mae$model %in% c(1:5),1:3] %>%  # *exclude M6 bc it's unfairly favoured in D5 and D8
  tidyr::pivot_wider(names_from = model, values_from = MAE)
wide_mae[, -1] <- t(apply(wide_mae[, -1], 1, function(x) {
  100 * (x - min(x)) / min(x) } ))
wide_mae <- round(as.data.frame(wide_mae), 2)  #percent change in MAE from minimum
colMeans(wide_mae[wide_mae$data %in% c(1,2,5,6),])  #when effort effect is proportional
colMeans(wide_mae[wide_mae$data %in% c(3,4,7,8),])  #when effort effect is threshold or constant

# MAE > high Effort values
wide_mae_max <- save_mae[save_mae$model %in% c(1:5),c(1,2,4)] %>%  # *exclude M6 bc it's unfairly favoured in D5 and D8
  tidyr::pivot_wider(names_from = model, values_from = MAE_max)
wide_mae_max[, -1] <- t(apply(wide_mae_max[, -1], 1, function(x) {
  100 * (x - min(x)) / min(x) } ))
wide_mae_max <- round(as.data.frame(wide_mae_max), 2)  #percent change in MAE from minimum
colMeans(wide_mae_max[wide_mae_max$data %in% c(1,2,5,6),])  #when proportional
colMeans(wide_mae_max[wide_mae_max$data %in% c(3,4,7,8),])  #when threshold or constant

# MAE < low Effort values
wide_mae_min <- save_mae[save_mae$model %in% c(1:5),c(1,2,5)] %>%  # *exclude M6 bc it's unfairly favoured in D5 and D8
  tidyr::pivot_wider(names_from = model, values_from = MAE_min)
wide_mae_min[, -1] <- t(apply(wide_mae_min[, -1], 1, function(x) {
  100 * (x - min(x)) / min(x) } ))
wide_mae_min <- round(as.data.frame(wide_mae_min), 2)  #percent change in MAE from minimum
colMeans(wide_mae_min[wide_mae_min$data %in% c(1,2,5,6),])  #when proportional
colMeans(wide_mae_min[wide_mae_min$data %in% c(3,4,7,8),])  #when threshold or constant

# Estimates - did they recover the true effect?
# These summaries are = 1 if 95% CI overlaps the true estimate, and = 0 otherwise
estimate_summary <- estimate_accuracy(results = results)
estimate_summary[estimate_summary$model==1,]  #examine each model

# Smoothers - did they recover the true effect?
prop_effort_smoothers()  #proportional and constant effort effects: red line should be within CIs
thresh_effort_smoothers()  #threshold effort effects: red line should be within CIs
temp_effort_smoothers()  #the temp smoother: red line should be within CIs; visually inspect domed effects
