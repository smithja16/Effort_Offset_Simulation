##############################################################
####   Code to generate and evaluate catch data...        ####
####   ... given Site, Temperature, and Effort effects... ####
####   ... to test how to best model Effort in GLMs.      ####
####   J.A.Smith NSW DPI 30/8/24                          ####
##############################################################

## This script generates an example figure showing numerous effort marginals
## The data generated are here are simpler than in the other function

# Load libraries
library(ggplot2)
library(mgcv)
library(gridExtra)

#set.seed(117)

# Generate data
generate_data <- function(n = 200, type = "proportional") {
  effort <- runif(n, 1, 10)
  
  if (type == "proportional") {
    true_lambda <- 2 * effort
    lambda <- true_lambda + rnorm(n, 0, 0.5)
  } else {
    true_lambda <- pmin(2 * effort, 10)
    lambda <- true_lambda + rnorm(n, 0, 0.5)
  }
  
  count <- rpois(n, lambda)
  #count <- rnbinom(n, mu = lambda, size = 1.1)
  data.frame(count = count, effort = effort, true_lambda = true_lambda)
}

# Generate datasets
data_prop <- generate_data(type = "proportional")
data_thresh <- generate_data(type = "threshold")

# Function to fit models and extract predictions
fit_models <- function(data) {
  models <- list(
    m1 = glm(count ~ offset(log(effort)), family = poisson, data = data),
    m2 = glm(count ~ log(effort), family = poisson, data = data),
    m3 = glm(count ~ effort, family = poisson, data = data),
    m4 = gam(count ~ s(log(effort)), family = poisson, data = data),
    m5 = gam(count ~ s(effort), family = poisson, data = data),
    m6 = glm(count ~ offset(effort), family = poisson, data = data)
  )
  new_data <- data.frame(effort = seq(min(data$effort), max(data$effort), length.out = 100))
  predictions <- lapply(models, function(m) {
    cbind(new_data, count = predict(m, newdata = new_data, type = "response"))
  })
  list(models = models, predictions = predictions)
}

# Fit models
fits_prop <- fit_models(data_prop)
fits_thresh <- fit_models(data_thresh)

# Create plots
create_plot <- function(data, fits, title, true_relationship) {
  ggplot(data, aes(x = effort, y = count)) +
    geom_point(alpha = 0.4, size = 0.7) +
    geom_line(aes(y = true_lambda, color = "True Relationship"), linewidth = 1) +
    geom_line(data = fits$predictions$m1, aes(color = "GLM: offset(log(effort))"), linetype = "solid") +
    geom_line(data = fits$predictions$m2, aes(color = "GLM: log(effort)"), linetype = "dashed") +
    geom_line(data = fits$predictions$m3, aes(color = "GLM: effort"), linetype = "dotdash") +
    geom_line(data = fits$predictions$m4, aes(color = "GAM: s(log(effort))"), linetype = "longdash") +
    geom_line(data = fits$predictions$m5, aes(color = "GAM: s(effort)"), linetype = "dashed") +
    geom_line(data = fits$predictions$m6, aes(color = "GAM: s(effort)"), linetype = "dashed") +
    scale_color_manual(values = c("True Relationship" = "black",
                                  "GLM: offset(log(effort))" = "red",
                                  "GLM: log(effort)" = "blue",
                                  "GLM: effort" = "darkgreen",
                                  "GAM: s(log(effort))" = "purple",
                                  "GAM: s(effort)" = "orange",
                                  "GLM: offset(effort)" = "pink")) +
    ylim(0, max(data)+2) +
    labs(title = title, x = "Effort", y = "Count", color = "Models") +
    theme_minimal() +
    theme(legend.position = "none")
}

plot_prop <- create_plot(data_prop, fits_prop, "Proportional Effort-Abundance", "2 * effort")
plot_thresh <- create_plot(data_thresh, fits_thresh, "Threshold Effort-Abundance", "min(2 * effort, 10)")

# Combine plots
combined_plot <- grid.arrange(plot_prop, plot_thresh, ncol = 2)

# Create a function to add a common legend
add_legend <- function(g1, g2, legend_title) {
  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  legend <- g_legend( ggplot() +
      geom_line(aes(x = 1, y = 1, color = "True Relationship"), linewidth = 1) +
      geom_line(aes(x = 1, y = 1, color = "GLM: offset(log(effort))"), linetype = "solid") +
      geom_line(aes(x = 1, y = 1, color = "GLM: log(effort)"), linetype = "dashed") +
      geom_line(aes(x = 1, y = 1, color = "GLM: effort"), linetype = "dotdash") +
      geom_line(aes(x = 1, y = 1, color = "GAM: s(log(effort))"), linetype = "longdash") +
      geom_line(aes(x = 1, y = 1, color = "GAM: s(effort)"), linetype = "dashed") +
      geom_line(aes(x = 1, y = 1, color = "GLM: offset(effort)"), linetype = "dashed") +
      scale_color_manual(values = c("True Relationship" = "black",
                                    "GLM: offset(log(effort))" = "red",
                                    "GLM: log(effort)" = "blue",
                                    "GLM: effort" = "green",
                                    "GAM: s(log(effort))" = "purple",
                                    "GAM: s(effort)" = "orange",
                                    "GLM: offset(effort)" = "pink")) +
      theme(legend.position = "bottom") +
      labs(color = legend_title) )
  
  grid.arrange(arrangeGrob(g1, g2, ncol = 2),
               legend,
               nrow = 2,
               heights = c(10, 1))
}

# Create the final plot with legend
final_plot <- add_legend(plot_prop, plot_thresh, "Models")
