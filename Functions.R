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


## Code to extract and plot some results
## I wrote this as a function to keep the main script cleaner

save_plot_results <- function(ylims) {
  
  # Write data frames to save data
  # Estimates and SEs (don't save smoothers here)
  save_coefs <- expand.grid(data = 1:8, model=1:6)
  save_coefs$Intercept <- NA; save_coefs$Temperature <- NA
  save_coefs$SiteB <- NA; save_coefs$SiteC <- NA
  save_coefs$logEffort <- NA; save_coefs$Effort <- NA
  save_SEs <- save_coefs
  
  # MAE (or could use RMSE)
  save_mae <- expand.grid(data = 1:8, model=1:6)
  save_mae$MAE <- NA  #all data
  save_mae$MAE_max <- NA  #high effort value data
  save_mae$MAE_min <- NA  #low effort value data
  
  par(mfrow=c(4,2), mar=c(2.5,3,2,1))
  colx <- c("red","blue","darkgreen","purple","green","orange")
  ylims <- ylims  #ylim for plot of each data scenario
  for (dd in 1:8) {  #data scenarios
    
    resultsx <- get(paste0("results_",dd))
    datax <- get(paste0("data",dd))
    
    for (mm in 1:6) {  #model types
      
      Mx <- resultsx[[paste0("M",mm)]]
      rowx <- which(save_mae$data==dd & save_mae$model==mm)
      
      # Estimates and SEs
      coefsx <- Mx$coefficients[1:5]
      SEsx <-  summary(Mx)$se[1:5]
      clean_names <- gsub("[()]", "", names(coefsx))
      clean_names <- clean_names[clean_names %in% names(save_coefs)]
      save_coefs[rowx,clean_names] <- coefsx[which(clean_names %in% names(save_coefs))]
      save_SEs[rowx,clean_names] <- SEsx[which(clean_names %in% names(save_SEs))]
      
      # MAE, in this case > 15 units is high effort, < 8 units low effort
      save_mae$MAE[rowx] <- MAE(Mx$fitted.values, datax$Count)
      save_mae$MAE_max[rowx] <- MAE(Mx$fitted.values[datax$Effort > 15], datax$Count[datax$Effort > 15])
      save_mae$MAE_min[rowx] <- MAE(Mx$fitted.values[datax$Effort < 8], datax$Count[datax$Effort < 8])
      
      # Marginal effort effects
      effort_seq <- seq(4, 20, length.out = 100)
      log_effort_seq <- log(effort_seq)
      new_data <- data.frame(Site="A",
                             Temperature=20,
                             Effort=effort_seq,
                             logEffort=log_effort_seq)
      pred <- predict(Mx, newdata=new_data, type="response")
      
      # Plot the marginals (this looks complicated due to necessary changes in how the 'true effect' is plotted)
      if (mm == 1) {
        plot(effort_seq, pred, type='l', ylim=c(0,ylims[dd]), xlim=c(4,20),
             xlab="Effort", ylab="Count", lwd=2, las=1,
             main=paste0("Effort marginal, data=",dd), col=colx[mm])
        legend("topleft", legend=c("model1","model2","model3","model4","model5","model6"),
               col=colx, lty=1:6, lwd=2, cex=0.8)
        if (dd %in% c(1,2,6)) {  #add proportional 'true effect'
          true_effect <- exp(0 + 0 + 0.1 * mean(datax$Temperature)) * effort_seq
          lines(effort_seq, true_effect, lwd=3, col="darkgrey")
        }
        if (dd == 1) { points(data1$Effort, data1$Count, col="grey", pch=20, cex=0.8) }
        if (dd %in% c(3,4,8)) {  #add threshold 'true effect'
          effort_threshold <- 12  #same as in the 'generate_counts' function
          median_effort <- median(datax$Effort)
          save_preds <- rep(0,6)
          for (mm in 1:6) {  #get mean prediction from all models
            Mxx <- resultsx[[paste0("M",mm)]]
            predx <- predict(Mxx, newdata=new_data, type="response")
            save_preds[mm] <- predx[which.min(abs(new_data$Effort - median_effort))]
          }
          rescale_factor <- mean(save_preds) / (min(median_effort, effort_threshold))
          true_effect <- rescale_factor * pmin(effort_seq, effort_threshold)
          lines(effort_seq, true_effect, lwd=3, col="darkgrey")
        }
        if (dd == 7) {
          true_effect <- mean(datax$Count)
          lines(effort_seq, rep(true_effect, length(effort_seq)), lwd=3, col="darkgrey")
        }
      } else {
        lines(effort_seq, pred, col=colx[mm], lty=mm, lwd=2)
      }
    }
  }
  return(list(save_coefs = save_coefs,
         save_SEs = save_SEs,
         save_mae = save_mae))
}


## Code to measure whether the 95% CI of an estimate encompasses the true value
## I wrote this as a function to keep the main script cleaner

estimate_accuracy <- function(results) {
  
  save_coefs <- results$save_coefs
  save_SEs <- results$save_SEs
  summary_in_out <- save_coefs
  for (vv in c("Intercept","Temperature","SiteB","SiteC","logEffort","Effort")) {
    if (vv == "Intercept") { target <- 0 }  #the true values
    if (vv == "Temperature") { target <- 0.1 }
    if (vv == "SiteB") { target <- 0.5 }
    if (vv == "SiteC") { target <- -0.5 }
    if (vv == "logEffort") { target <- 1 }
    if (vv == "Effort") { target <- 0 }
    for (rr in 1:nrow(save_coefs)) {
      coefx <- save_coefs[rr,vv]
      SEx <- save_SEs[rr,vv]
      if (is.na(coefx) == F) {
        lower_ci <- coefx - SEx*1.96
        upper_ci <- coefx + SEx*1.96
        if (lower_ci <= target && upper_ci >= target) {
          summary_in_out[rr,vv] <- 1
        } else {
          summary_in_out[rr,vv] <- 0
        }
      }
    }
  }
  return(summary_in_out)
  
}


## Code to plot and evaluate whether the models with effort smoothers...
## ...accurately recover the proportional effort effects
## I wrote this as a function to keep the main script cleaner

prop_effort_smoothers <- function() {
  
  # Evaluate Effort smoother in M4, M5 & M6
  # Successful if red line (prop. Effort effect is within CIs of GAM)
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  for (dd in c(1,2,5,6,7)) {
    resultsx <- get(paste0("results_",dd))
    for (mm in c(4,6)) {
      Mx <- resultsx[[paste0("M",mm)]]
      px <- ggeffects::predict_response(Mx, terms = c("logEffort"), 
                                        condition = c(Temperature=20, Site="A"))
      plot(exp(px$x), px$predicted, type="l",
           ylim=c(0,max(px$conf.high)), main=paste0("M",mm," - Data ",dd), las=1,
           ylab="Count", xlab="Effort")
      lines(exp(px$x), px$conf.low, lty=2)
      lines(exp(px$x), px$conf.high, lty=2)
      if (dd == 7) { #no effort effect - slope = 0
        idx <- which.min(abs(px$x - log(12)))  #find a point near middle that GAM response runs through
        mid_pointx <- px$x[idx]
        mid_pointy <- px$predicted[idx]
        lines(exp(px$x), rep(mid_pointy, length(px$x)), col="red")
      } else {
        idx <- which.min(abs(px$x - log(12)))  #find a point near middle that GAM response runs through
        mid_pointx <- px$x[idx]
        mid_pointy <- px$predicted[idx]
        intercept <- log(mid_pointy) - log(12) - 0.1*20
        lines(exp(px$x), exp(intercept + 0.1*20 + px$x), col="red")
      }
    }
    for (mm in c(5)) {
      Mx <- resultsx[[paste0("M",mm)]]
      px <- ggeffects::predict_response(Mx, terms = c("Effort"), 
                                        condition = c(Temperature=20, Site="A"))
      plot(px$x, px$predicted, type="l",
           ylim=c(0,max(px$conf.high)), main=paste0("M",mm," - Data ",dd), las=1,
           ylab="Count", xlab="Effort")
      lines(px$x, px$conf.low, lty=2)
      lines(px$x, px$conf.high, lty=2)
      if (dd == 7) { #no effort effect - slope = 0
        idx <- which.min(abs(px$x - 12))  #find a point near middle that GAM response runs through
        mid_pointx <- px$x[idx]
        mid_pointy <- px$predicted[idx]
        lines(px$x, rep(mid_pointy, length(px$x)), col="red")
      } else {
        idx <- which.min(abs(px$x - 12))  #find a point near middle that GAM response runs through
        mid_pointx <- px$x[idx]
        mid_pointy <- px$predicted[idx]
        intercept <- log(mid_pointy) - log(12) - 0.1*20
        lines(px$x, exp(intercept + 0.1*20 + log(px$x)), col="red")
      }
    }
  }
  
}


## Code to plot and evaluate whether the models with effort smoothers...
## ...accurately recover the threshold effort effects
## I wrote this as a function to keep the main script cleaner

thresh_effort_smoothers <- function() {
  
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  for (dd in c(3,4,8)) {
    datax <- get(paste0("data",dd))
    resultsx <- get(paste0("results_",dd))
    for (mm in c(4,6)) {
      Mx <- resultsx[[paste0("M",mm)]]
      px <- ggeffects::predict_response(Mx, terms = c("logEffort"), 
                                        condition = c(Temperature=20, Site="A"))
      plot(exp(px$x), px$predicted, type="l",
           ylim=c(0,max(px$conf.high)), main=paste0("M",mm," - Data ",dd), las=1,
           ylab="Count", xlab="Effort")
      lines(exp(px$x), px$conf.low, lty=2)
      lines(exp(px$x), px$conf.high, lty=2)
      effort_threshold <- 12  #previously max(datax$Effort)/2
      true_effect <- exp(0 + 0.1 * 20) * pmin(exp(px$x), effort_threshold)
      if (dd == 8) {  #needs rescaling due to strong collinearity
        median_effort <- median(datax$Effort)
        median_pred <- px$predicted[which.min(abs(exp(px$x) - median_effort))]
        rescale_factor <- median_pred / (min(median_effort, effort_threshold))
        true_effect <- rescale_factor * pmin(exp(px$x), effort_threshold)
      }
      lines(exp(px$x), true_effect, col="red")
    }
    for (mm in c(5)) {
      Mx <- resultsx[[paste0("M",mm)]]
      px <- ggeffects::predict_response(Mx, terms = c("Effort"), 
                                        condition = c(Temperature=20, Site="A"))
      plot(px$x, px$predicted, type="l",
           ylim=c(0,max(px$conf.high)), main=paste0("M",mm," - Data ",dd), las=1,
           ylab="Count", xlab="Effort")
      lines(px$x, px$conf.low, lty=2)
      lines(px$x, px$conf.high, lty=2)
      effort_threshold <- 12  #previousy max(datax$Effort)/2
      true_effect <- exp(0 + 0.1 * 20) * pmin(px$x, effort_threshold)
      if (dd == 8) {
        median_effort <- median(datax$Effort)
        median_pred <- px$predicted[which.min(abs(px$x - median_effort))]
        rescale_factor <- median_pred / (min(median_effort, effort_threshold))
        true_effect <- rescale_factor * pmin(px$x, effort_threshold)
      }
      lines(px$x, true_effect, col="red")
    }
  }
  
}


## Code to plot and evaluate whether the model with a temperature smoother...
## ...could accurately recover the temperature effect
## I wrote this as a function to keep the main script cleaner

temp_effort_smoothers <- function() {
  
  par(mfrow=c(2,2))
  for (dd in c(1:4,6:7)) {
    resultsx <- get(paste0("results_",dd))
    Mx <- resultsx[["M6"]]
    px <- ggeffects::predict_response(Mx, terms = c("Temperature"), 
                                      condition = c(logEffort=log(12), Site="A"))
    plot(px$x, px$predicted, type="l", ylim=c(0,max(px$conf.high)),
         main=paste("M6 - Data",dd), las=1, ylab="Count", xlab="Temperature")
    lines(px$x, px$conf.low, lty=2)
    lines(px$x, px$conf.high, lty=2)
    idx <- which.min(abs(px$x - 20))  #find a point near middle that GAM response runs through
    mid_pointx <- px$x[idx]
    mid_pointy <- px$predicted[idx]
    intercept <- log(mid_pointy) - log(12) - 0.1*mid_pointx
    lines(px$x, exp(intercept + 0.1*px$x + log(12)), col="red")
  }
  for (dd in c(5,8)) {  #visually inspect Temp smoothers
    resultsx <- get(paste0("results_",dd))
    Mx <- resultsx[["M6"]]
    plot(Mx, select=1, main=paste0("M6, Data",dd))
  }
  
}
