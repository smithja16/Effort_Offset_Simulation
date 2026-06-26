##############################################################
####   Code to generate and evaluate catch data...        ####
####   ... given Site, Temperature, and Effort effects... ####
####   ... to test how to best model Effort in GLMs.      ####
####   J.A.Smith NSW DPI - June 2026                      ####
##############################################################

## This script generates an example figure showing how when effort
## is endogenous to abundance, this can confound the effort-abundance
## relationship, and bias estimated trends in abundance. This is 
## put in a fisheries context due to the suitability of this context 
## and importance of this issue for catch rate standardizations.

## Scenario:
## The Effort-Abundance relationship is proportional, but also:
## 1) When fish abundance is higher, catch rates are higher...
## 2) When catch rates are higher, fishers respond by increasing effort...
## 3) This creates a positive correlation between effort and abundance. 
## This means that in the years where abundance is higher, and catch
## rate is higher, we also have higher effort. 

## The outcome of this scenario is that the offset model tracks the
## true abundance trend more closely than the covariate model, because
## the latter attributes some of the abundance trend to the effort term.
## The covariate model "explains away" some of the abundance signal
## through the effort term, which flattens the abundance trend that
## doesn't show the full extent of the decline.

## This simulation result clearly supports the use of offsets
## rather than covariates when:
## The effort-catch relationship is truly proportional (or close to it) &
## Effort allocation responds to catch rates (targeting behavior)


# Load libraries
library(ggplot2)
library(reshape2)
library(patchwork)

set.seed(333)

# Parameters
n_years <- 10
n_locations <- 30  # 'Abundance' is measured at numerous locations
n_observations <- n_years * n_locations

# Create year and location variables
year <- rep(1:n_years, each = n_locations)
location <- rep(1:n_locations, times = n_years)

# True abundance trend (declining)
true_year_effect <- exp(-0.1 * (1:n_years))

# Location effects (random variation by location)
location_effect <- rlnorm(n_locations, meanlog = 0, sdlog = 0.2)

# True abundance for each observation
baseline_abundance <- 10  # Added this factor to increase overall abundance
true_abundance <- baseline_abundance * true_year_effect[year] * location_effect[location]

# Add small random noise to true abundance
abundance_with_noise <- true_abundance * rlnorm(n_observations, meanlog = 0, sdlog = 0.1)

# Catchability coefficient
q <- 0.5  #high enough to avoid too many zeros

# Expected CPUE if effort were constant
expected_cpue <- q * abundance_with_noise

# The feedback loop:
# 1. Base effort (before feedback) - moderate variation
base_effort <- exp(rnorm(n_observations, mean = 0, sd = 0.2))

# 2. Adjust effort based on expected CPUE - strong feedback
# Higher CPUE -> higher effort (clear positive relationship)
effort_response <- function(cpue_values) {
  # create a strong positive relationship between CPUE and effort
  normalized_cpue <- (cpue_values - min(cpue_values)) / 
    (max(cpue_values) - min(cpue_values))
  # Convert to effort multiplier (range ~0.5 to 3)
  multiplier <- 0.5 + 2.5 * normalized_cpue
  return(multiplier)
}

effort_multiplier <- effort_response(expected_cpue)
effort <- base_effort * effort_multiplier

# Generate catch with this effort
expected_catch <- q * abundance_with_noise * effort
catch <- rpois(n_observations, expected_catch)

# Calculate observed CPUE
cpue <- catch / effort

data <- data.frame(
  year = factor(year),
  location = factor(location),
  effort = effort,
  catch = catch,
  cpue = cpue,
  true_abundance = abundance_with_noise,
  true_year_effect = true_year_effect[year] )

# Check for zeros and summarize catches by year
zero_catch_count <- sum(data$catch == 0)
print(paste("Number of zero catches:", zero_catch_count, "out of", n_observations, 
            "observations (", round(100*zero_catch_count/n_observations, 1), "%)"))

catch_by_year <- aggregate(data$catch ~ data$year, FUN=sum)
print("Catch totals by year:")
print(catch_by_year)

# Verify the relationship between abundance and effort
p1 <- ggplot(data, aes(x = true_abundance, y = effort)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(title = "a) True Abundance vs Effort",
       x = "True Abundance", 
       y = "Effort")
print(p1)

# Plot the specified time series of abundance
mean_n <- aggregate(true_abundance ~ year, FUN = mean)
p2 <- ggplot(mean_n, aes(x = year, y = true_abundance)) +
  geom_line(linewidth = 1, col = "#619CFF") + 
  labs(title = "b) True Abundance decline",
       x = "Year", 
       y = "True Abundance")
print(p2)


# Fit catch rate standardization models
# 1. Effort as offset
m_offset <- glm(catch ~ year + offset(log(effort)), family = quasipoisson, data = data)

# 2. Effort as covariate
m_covariate <- glm(catch ~ year + log(effort), family = quasipoisson, data = data)

# Extract year coefficients
get_abundance_index <- function(model, reference_year = 1) {
  year_coefs <- coef(model)[grep("year", names(coef(model)))]
  # Add reference year (which has coefficient 0)
  year_coefs <- c(0, year_coefs)
  # Convert to abundance index
  abundance_index <- exp(year_coefs)
  # Scale to first year = 1
  abundance_index <- abundance_index / abundance_index[reference_year]
  return(abundance_index)
}

# Get abundance indices
index_offset <- get_abundance_index(m_offset)
index_covariate <- get_abundance_index(m_covariate)

# True index (normalized to year 1)
true_index <- true_year_effect / true_year_effect[1]

# Calculate error metrics
rmse_offset <- sqrt(mean((index_offset - true_index)^2))
rmse_covariate <- sqrt(mean((index_covariate - true_index)^2))

# Print effort coefficient from covariate model
print(paste("Effort coef. in covariate model:", coef(m_covariate)["log(effort)"]))

# Print error metrics
print(paste("RMSE for offset model:", rmse_offset))
print(paste("RMSE for covariate model:", rmse_covariate))

# Plot
years <- 1:n_years
results_df <- data.frame(
  Year = rep(years, 3),
  Model = factor(rep(c("True", "Offset", "Covariate"), each = n_years)),
  Index = c(true_index, index_offset, index_covariate) )

p3 <- ggplot(results_df, aes(x = Year, y = Index, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "d) Comparison of Abundance Indices",
       subtitle = paste("RMSE: Offset:", round(rmse_offset, 4), 
                        ", Covariate:", round(rmse_covariate, 4)),
       y = "Relative Abundance Index",
       x = "Year") + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.8))
print(p3)



## The estimated effort coefficient from the covariate model
effort_coef <- coef(m_covariate)["log(effort)"]

# Sequence of effort values
effort_seq <- seq(min(data$effort), max(data$effort), length.out = 100)

# Set a constant abundance and show the functional form
# of how catch relates to effort under each model

# Calculate expected catch-effort relationships
# For the offset model, catch ~ effort^1
# For the covariate model, catch ~ effort^coefficient
df_relationships <- data.frame(
  effort = effort_seq,
  # Offset model (proportional relationship)
  Offset = effort_seq,
  # Covariate model (estimated relationship)
  Covariate = effort_seq^effort_coef )

# Scale both to start at the same point (for easier comparison)
df_relationships$Offset <- df_relationships$Offset / df_relationships$Offset[1]
df_relationships$Covariate <- df_relationships$Covariate / df_relationships$Covariate[1]

# Reshape and plot
df_long <- melt(df_relationships, id.vars = "effort", 
                variable.name = "relationship", value.name = "relative_catch")

p4 <- ggplot(df_long, aes(x = effort, y = relative_catch, color = relationship)) +
  geom_line(linewidth = 1.5) +
  labs(title = "c) Effort-Catch Relationships",
       subtitle = paste("Covariate model: effort coef =", round(effort_coef, 3),
                        "\nOffset model: effort coef = 1 (fixed)"),
       x = "Effort", 
       y = "Relative Expected Catch",
       color = "Model Type") +
  scale_color_manual(values = c("Offset" = "#00BA38", "Covariate" = "#F8766D")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.2,0.8))
print(p4)


(p1 + p2) / (p4 + p3)

