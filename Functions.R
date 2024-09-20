##############################################################
####   Functions to generate and evaluate catch data...   ####
####   ... given Site, Temperature, and Effort effects    ####
####   J.A.Smith NSW DPI 30/8/24                          ####
##############################################################

## Generic RMSE function
RMSE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  sqrt(mean((p - o)^2))
}


## Generic MAE function (less sensitive to outliers, so better choice in some cases)
MAE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  sum(abs(p - o))/length(p)
}


## *Workhorse data generation function*
## This generates count data from a negative binomial distribution with specified theta
## Counts are determined by specified effort, site (3 levels), and temperature effects
## The effort effect can be "proportional", "threshold", or "constant"...
## ... "proportional" is a 1:1 relationship assumed by an offset term
## ... "threshold" is a simple rising then constant threshold effect
## ... "constant" is no effort effect (count is independent of effort)
## The temp effect can be "linear" (uses beta_temp) or "nonlinear" (uses temp_optimum and temp_breadth)
## The Site and Temperature effects are arbitrary but scaled to be meaningful in magnitude
## All 3 variables can be correlated by specifying covariance in the MVN function
## This does mean the distribution of effort values is normal, which may not always be the case
## Adding correlation can change the magnitude of marginal effects and thus magnitude of counts

generate_counts <- function(n_obs = 1000,  #number of counts to generate 
                            temp_effort_cor = 0,  # covariance between Temp and Effort 
                            site_effort_cor = 0,  # covariance between Site and Effort
                            effort_type = c("proportional", "threshold", "constant"),
                            temp_effect = c("linear", "nonlinear"),
                            temp_optimum = 20,  # temp at which nonlinear effect peaks
                            temp_breadth = 5,   # controls how quickly nonlinear effect drops off
                            beta_temp = 0.1,    # coefficient for linear temp effect
                            seed = NULL) {  #can set seed for reproducibility
  
  set.seed(seed)
  
  effort_type <- match.arg(effort_type)
  temp_effect <- match.arg(temp_effect)
  
  ## Generate site variable
  site <- as.factor(sample(c("A", "B", "C"), n_obs, replace = TRUE))
  
  ## Generate temperature and effort, considering both correlations simultaneously
  if (temp_effort_cor == 0 && site_effort_cor == 0) {
    temperature <- rnorm(n_obs, mean = 20, sd = 3)
    effort <- rnorm(n_obs, mean = 12, sd = 2.5)
  } else {
    # Create correlation matrix
    cor_matrix <- matrix(c(1, temp_effort_cor, 0,
                           temp_effort_cor, 1, site_effort_cor,
                           0, site_effort_cor, 1), nrow = 3)
    
    # Generate multivariate normal distribution
    mvn_vars <- MASS::mvrnorm(n = n_obs, mu = c(0, 0, 0), Sigma = cor_matrix)
    
    # Extract and scale variables
    temperature <- mvn_vars[,1] * 3 + 20  # scale to desired sd = 3 and mean = 20
    effort_base <- mvn_vars[,2] * 2.5 + 12  # scale to desired sd = 2.5 and mean = 12
    site_effect <- mvn_vars[,3] * site_effort_cor * 2.3  # scale effect size (amplify the Site effect)
    
    # Combine base effort and site effect
    effort <- effort_base + site_effect
    
    # Adjust site effect based on actual site
    site_adjustments <- c(A = -mean(site_effect), 
                          B = site_effort_cor * 2.3 - mean(site_effect), 
                          C = -site_effort_cor * 2.3 - mean(site_effect))
    effort <- effort + site_adjustments[site]
  }
  
  ## Calculate linear predictor
  beta_site <- c(0, 0.5, -0.5)  # Effect sizes for sites B and C relative to A
  
  # Temperature effect based on user choice
  if (temp_effect == "linear") {
    temp_effect <- beta_temp * temperature
  } else {  # nonlinear
    temp_effect <- exp(-((temperature - temp_optimum)^2) / (2 * temp_breadth^2))
  }
  
  linear_pred <- beta_site[as.numeric(site)] + temp_effect
  
  ## Apply effort effect based on type
  if (effort_type == "proportional") {
    linear_pred <- linear_pred + log(effort)
  } else if (effort_type == "threshold") {
    effort_threshold <- 12  # previously used max(effort)/2, but this varied a lot due to random large values
    linear_pred <- linear_pred + log(pmin(effort, effort_threshold))
  } else if (effort_type == "constant") {
    # No effect of effort
  }
  
  ## Generate counts using negative binomial distribution or Poisson
  theta <- 3  # controls overdispersion; higher means less dispersion
  mu <- exp(linear_pred)
  counts <- rnbinom(n_obs, mu = mu, size = theta)
  #counts <- rpois(n_obs, lambda = mu)  #useful test for increasing signal:noise
  
  # Create and return a data frame
  return(data.frame(Count = counts,
                    Site = site,
                    Temperature = temperature,
                    Effort = effort))
}
