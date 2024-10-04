##############################################################
####   Code to generate and evaluate catch data...        ####
####   ... given Site, Temperature, and Effort effects... ####
####   ... to test how to best model Effort in GLMs.      ####
####   J.A.Smith NSW DPI 30/8/24                          ####
##############################################################

## This script generates an example figure showing the difference in 
## assumed effort-probability relationship when adding an effort
## covariate to a binomial GLM (Figure S2)

## It generates biniomial data from a logistic driven by temperature and effort
## Then explore the relative shape of transformed and untransformed effort covariates
## The goal is not to show the most suitable for these data, just their relative shapes
## on link and original scales
## Here "capture" means sampling/detecting/observing etc


# Load libraries
library(ggplot2)
library(ggeffects)
library(patchwork)

#set.seed(117)

# Create fake data
n <- 1000
temperature <- rnorm(n, mean=20, sd=5)
effort <- runif(n, min=1, max =10)
logit_p <- -2 + 0.1*temperature + 0.3*effort
p <- 1 / (1 + exp(-logit_p))
capture <- rbinom(n, size=1, prob=p)

data <- data.frame(capture=capture, temperature=temperature, effort=effort)
data$logEffort <- log(data$effort)

# Fit binomial glms
model1 <- glm(capture ~ temperature + effort, family=binomial, data=data)
model2 <- glm(capture ~ temperature + logEffort, family=binomial, data=data)

#summary(model1)
#summary(model2)

# Calculate marginal effects
marg_temp1 <- predict_response(model1, terms="temperature")
marg_temp2 <- predict_response(model2, terms="temperature")
newdata <- data.frame(effort=seq(1,10,by=0.1), temperature=20.08)
marg_effort1 <- predict_response(model1, terms=newdata)
newdata <- data.frame(logEffort=log(seq(1,10, by=0.1)), temperature=20.08)
marg_effort2 <- predict_response(model2, terms=newdata)
marg_effort2$x <- exp(marg_effort2$x)  #transform back to original scale

# Combine
marg_effort1$model <- "Untransformed Effort"
marg_effort2$model <- "Log-transformed Effort"
marg_effort_combined <- rbind(marg_effort1, marg_effort2)

# Calculate on link scale
marg_effort1$predicted_link <- qlogis(marg_effort1$predicted)
marg_effort1$conf.low_link <- qlogis(marg_effort1$conf.low)
marg_effort1$conf.high_link <- qlogis(marg_effort1$conf.high)

marg_effort2$predicted_link <- qlogis(marg_effort2$predicted)
marg_effort2$conf.low_link <- qlogis(marg_effort2$conf.low)
marg_effort2$conf.high_link <- qlogis(marg_effort2$conf.high)

marg_effort_combined_link <- rbind(marg_effort1, marg_effort2)

# Plot marginal effects on original scale
plot_temp <- ggplot() +
  geom_line(data = marg_temp1, aes(x = x, y = predicted, color = "Untransformed Effort")) +
  geom_ribbon(data = marg_temp1, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = marg_temp2, aes(x = x, y = predicted, color = "Log-transformed Effort")) +
  geom_ribbon(data = marg_temp2, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(title = "Marginal Effect of Temperature",
       x = "Temperature",
       y = "Predicted Probability of Capture",
       color = "Model") +
  theme_minimal()

plot_effort <- ggplot(marg_effort_combined, aes(x = x, y = predicted, color = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model), alpha = 0.1) +
  labs(title = "Marginal Effect of Effort - original scale",
       x = "Effort",
       y = "Predicted Probability of Capture",
       color = "Model",
       fill = "Model") +
  theme_classic()

# Plot effort marginal effect on link scale
plot_effort_link <- ggplot(marg_effort_combined_link, aes(x = x, y = predicted_link, color = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low_link, ymax = conf.high_link, fill = model), alpha = 0.1) +
  labs(title = "Marginal Effect of Effort - link scale",
       x = "Effort",
       y = "Predicted Log-odds of Capture",
       color = "Model",
       fill = "Model") +
  theme_classic()

# Display plots
print(plot_temp)
print(plot_effort)
print(plot_effort_link)

plot_effort / plot_effort_link
