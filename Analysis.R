##############################################################
####   Code to generate and evaluate catch data...        ####
####   ... given Site, Temperature, and Effort effects... ####
####   ... to test how to best model Effort in GLMs.      ####
####   J.A.Smith NSW DPI 30/8/24                          ####
##############################################################


## Load libraries
library(mgcv) 

## Load functions
source("Functions.R")


## Generate the data scenarios
## Here we generate 8 based on useful changes in collinearity, effort effects, and temp effects

# Number of data points
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

# Plot histograms for data scenarios; note that collinearity can change the magnitude of counts
par(mfrow=c(3,3))
for (dd in 1:8) {
  datax <- get(paste0("data",dd))
  hist(datax$Count, breaks=30,
       main=paste0("Data=",dd))
}


## For each data type, fit these six models
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

for (dd in 1:8) {  #each data scenario
  assign(paste0("results_",dd), fit_models(data_type = dd))
}

